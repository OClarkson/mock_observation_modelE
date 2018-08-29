#=============================================================================
# Module
#=============================================================================
import numpy  as np
import sys
from copy import deepcopy

import cgs
import ray_trace
import projection
import opacity
from molecules import *
import util_interp
import util_errors

TAU_MAX = 1000.


#=============================================================================
def raytrace_opacity( grid_wn, theta, d_theta, dict_atmprof, dict_griddata_logXSofWNTP, dict_NonCondensableGas, z_top, z_num, b_num, d_l, dict_geom, flag_molabs=True, flag_rayleigh=True, flag_cld=False ) :

    #------------------------------------------------
    # set up atmospheric profile
    #------------------------------------------------

    # interpolate atmospheric profile
    dict_atmprof_funcZ = {}
    param_x = deepcopy( dict_atmprof['z'] )
    for key in dict_atmprof :
        if key != 'z' :
            param_y = dict_atmprof[key]
            if ( np.any( param_y <= 0. ) or ( np.max( param_y ) / np.min( param_y ) ) < 10. ) :
                dict_atmprof_funcZ[key] = util_interp.interp_1d_boundary( param_x, param_y, logx=False, logy=False, ext=3 )
            else :
                dict_atmprof_funcZ[key] = util_interp.interp_1d_boundary( param_x, param_y, logx=False, logy=True, ext=3 )

                
    # initial set-up
    layer_b = np.linspace( 0, z_top, b_num )
    layer_z = np.linspace( 0., z_top, z_num ) # 1 km -> z_top
    # layer_z         = np.logspace( np.log10(1e5), np.log10(z_top), z_num ) # 1 km -> z_top

    b_bottom_indx, matrixB_z_min, matrixB_deflect, listB_func_lofz = ray_trace.refraction( layer_b, dict_geom['r_planet'], dict_atmprof_funcZ, d_l )

    # check the maximum pressure
    if len( dict_griddata_logXSofWNTP ) > 0 : 
        max_pres_probed = dict_atmprof_funcZ['plm']( matrixB_z_min[b_bottom_indx] )
        max_pres_in_xstbl = np.exp( np.max( dict_griddata_logXSofWNTP['coords'][:,1] ) )
        print 'max_pres_probed', max_pres_probed/cgs.mbar_to_barye
        print 'max_pres_in_xstbl', max_pres_in_xstbl/cgs.mbar_to_barye
        if ( max_pres_probed > max_pres_in_xstbl ) :
            util_errors.warning_longmsg( [ 'Maximum atmospheric pressure probed is larger than maximum pressure in lookuptable.' , 
                                           'Pressure higher than ' + str( max_pres_in_xstbl/cgs.mbar_to_barye ) + ' mbar is ignored.' ] )

    matrixB_area = np.zeros( b_num )
    matrixB_area[b_bottom_indx:] = projection.area_on_stellardisk( theta, d_theta, layer_b[b_bottom_indx:], matrixB_deflect[b_bottom_indx:], dict_geom ) 
#    matrixB_area[b_bottom_indx:] = projection.area_on_stellardisk( theta, d_theta, layer_b[b_bottom_indx:], matrixB_deflect[b_bottom_indx:], dict_geom )

    # opacity
    matrixWZ_nXS = np.zeros( [ len( grid_wn ), z_num ] )
    if flag_molabs and ( len( dict_griddata_logXSofWNTP ) > 0 ) :
        matrixWZ_nXS += opacity.get_nXS_molabs( layer_z, grid_wn, dict_griddata_logXSofWNTP, dict_NonCondensableGas, dict_atmprof_funcZ )
    if flag_rayleigh :
        matrixWZ_nXS += opacity.get_nXS_Rayleigh( layer_z, grid_wn, dict_NonCondensableGas, dict_atmprof_funcZ )
    if flag_cld :
        matrixWZ_nXS += opacity.get_nXS_cld( layer_z, grid_wn, dict_atmprof_funcZ )
    # delta l
    matrixBZ_dl      = np.zeros( [ b_num, z_num ] )

    for ii in xrange( b_bottom_indx, b_num - 1 ):

        marker_z = layer_z[np.where( layer_z > matrixB_z_min[ii] )]
        marker_l = listB_func_lofz[ii]( marker_z )
        marker_l_edge = ( marker_l[1:] + marker_l[:-1] ) / 2.
        marker_l_edge = np.r_[ 0., marker_l_edge, marker_l[-1] ]
        marker_deltal = np.diff( marker_l_edge )
        matrixBZ_dl[ii][np.where( layer_z > matrixB_z_min[ii] )[0][0]:] = marker_deltal

    matrixWB_tau = np.dot( matrixWZ_nXS, matrixBZ_dl.T )
    matrixWB_tau = np.clip( matrixWB_tau, 0, TAU_MAX )

#    matrixW_Ftransmit = np.dot( np.exp( -2. * matrixWB_tau ), matrixB_area )
    matrixW_Ftransmit = np.dot( np.exp( - 2. * matrixWB_tau ), matrixB_area )

    return matrixW_Ftransmit


