import netCDF4
import numpy as np
from copy import deepcopy

from molecules import *
import set_O3

from setup import l_O3, s_O3file

import cgs
import constants
import util_errors
import sys

#=============================================================================
def extract_limbprof( infile, dict_NonCondensableGas, p_min, p_max, g_planet ) :

    print "Reading " + infile + " for atmospheric profile..."

    # molecular weight of dry air
    mu_air_dry = 0.
    sum        = 0.
    weight_otherwise = 0.
    for molename in dict_NonCondensableGas :
        if dict_NonCondensableGas[molename]=='otherwise' :
            weight_otherwise = molecules[molename]['weight']
        else :
            mu_air_dry += dict_NonCondensableGas[molename]*molecules[molename]['weight']
            sum        += dict_NonCondensableGas[molename]
    mu_air_dry += ( 1. - sum )*weight_otherwise

    # read netCDF file
    ncfile_r = netCDF4.Dataset( infile, 'r', format='NETCDF3_64BIT' )

    plm = ncfile_r.variables['plm']
    lat = ncfile_r.variables['lat']
    lon = ncfile_r.variables['lon']

    n_lat = len( lat ) - 1
    n_lon = len( lon )


    # check the location of sub-stellar point
    tsurf = ncfile_r.variables['TempL'][0,:,:]
    lat_mesh, lon_mesh = np.meshgrid( lat, lon, indexing='ij' )
    lat_mesh_flat = lat_mesh.flatten()
    lon_mesh_flat = lon_mesh.flatten()
    lat_maxtsurf = lat_mesh_flat[np.argmax( tsurf )]
    lon_maxtsurf = lon_mesh_flat[np.argmax( tsurf )]
    if ( not np.fabs( lat_maxtsurf ) < 30. ) or ( not ( np.fabs( lon_maxtsurf ) > 150. ) ) :
        util_errors.warning_longmsg( [ 'This program assumes that the sub-stellar point is at 0 degree  ', 
                                       'in latitude and at -180 degree in longitude.                    ', 
                                       'Make sure your aijl file is consistent with it.                 ', 
                                       'Temperature maximum is found away from it hence this warning.   ', 
                                       'Ignore this warning if the planet is NOT synchronously rotating.'] )

    # west limb ( west from the substellar point,  90 degree longitude )
    ilon_limb_w = np.array( [ np.where(lon[:] > 90.0)[0][0] , np.where(lon[:] > 90.0)[0][0]-1 ] )
    lon_limb_w  = np.array( [ lon[ilon_limb_w[0]] , lon[ilon_limb_w[1]] ] )
    weight_w    = np.array( [ (lon_limb_w[1]-90.0)/(lon_limb_w[1]-lon_limb_w[0]) , (90.0-lon_limb_w[0])/(lon_limb_w[1]-lon_limb_w[0]) ] )

    # east limb ( east from the substellar point, -90 degree longitude )
    ilon_limb_e = np.array( [ np.where(lon[:] > -90.0)[0][0] , np.where(lon[:] > -90.0)[0][0]-1 ] )
    lon_limb_e  = np.array( [ lon[ilon_limb_e[0]] , lon[ilon_limb_e[1]] ] )
    weight_e    = np.array( [ (lon_limb_e[1]+90.0)/(lon_limb_e[1]-lon_limb_e[0]) , (-90.0-lon_limb_e[0])/(lon_limb_e[1]-lon_limb_e[0]) ] )
        
    ilon_limb_ew  = [ilon_limb_e ,  ilon_limb_w ]
    lon_limb_ew   = [ lon_limb_e,    lon_limb_w ]
    weight_ew     = [   weight_e,      weight_w ]
    latsign_ew    = [         1.,            -1.]
    grid_ew       = [ n_lat - np.arange( n_lat ), np.arange( n_lat ) ]

    #--------------------------------------------------------------------------
    # parameters to read
    param  = { 'z'      :{'scale': cgs.m_to_cm, 'offset':   0.   , 'fill' : 0. }, 
               'TempL'  :{'scale': 1e0        , 'offset':   0. , 'fill' : 0. }, 
               'q'      :{'scale': 1e0        , 'offset':   0.   , 'fill' : 0. }, 
               'cf'     :{'scale': 1e-2       , 'offset':   0.   , 'fill' : 0. }, 
               'icecld' :{'scale': 1e0        , 'offset':   0.   , 'fill' : 0. }, 
               'wtrcld' :{'scale': 1e0        , 'offset':   0.   , 'fill' : 0. } }
    #--------------------------------------------------------------------------

    list_dict_atmprof = []
    list_theta        = []

    # east limb or west limb
    for i_ew in xrange(2):

        # at different latitude
        for i_lat in grid_ew[i_ew] :

            if i_ew == 0 :
                theta = ( 90. - lat[i_lat] ) / 180. * np.pi
            else :
                theta = 2. * np.pi - ( 90. - lat[i_lat] ) / 180. * np.pi
            list_theta.append( theta )

            dict_atmprof = {}

            # input parameters
            for key in param:

                # masking
                filled_data = np.ma.filled( ncfile_r[key][:,:,:].astype('float'), param[key]['fill'] )

                dict_atmprof[key] =   weight_ew[i_ew][0] * filled_data[:,i_lat,ilon_limb_ew[i_ew][0]] \
                                    + weight_ew[i_ew][1] * filled_data[:,i_lat,ilon_limb_ew[i_ew][1]]
                
                # offset
                dict_atmprof[key] = dict_atmprof[key]*param[key]['scale'] + param[key]['offset']

            # unit conversion
            dict_atmprof['plm']   = plm[:]*cgs.mbar_to_barye
            # if ( dict_atmprof['plm'][0] > p_max ): 
            #     util_errors.warning_longmsg( [ 'Maximum pressure of atmosphere is larger than minimum pressure in lookuptable.' , 
            #                                    'Pressure higher than ' + str(p_max/cgs.mbar_to_barye) + 'bar is ignored.' ] )
            #     dict_atmprof['plm'][np.where( dict_atmprof['plm'] > p_max )] = p_max

            # extrapolate or truncate to p_min
            if ( dict_atmprof['plm'][-1] > p_min ):
                dict_atmprof = extrapolate_to_p_min( p_min, dict_atmprof, g_planet, param, mu_air_dry )
            else : 
                dict_atmprof = truncate_to_p_min( p_min, dict_atmprof, param )

            # input O3
            if l_O3 :
                dict_atmprof['xO3']   = set_O3.read_O3file( s_O3file, dict_atmprof['plm']/cgs.mbar_to_barye )

            # H2O vapor mixing ratio
            dict_atmprof['xH2O'] = dict_atmprof['q'] / constants.MU_H2O / ( dict_atmprof['q'] / constants.MU_H2O + ( 1 - dict_atmprof['q'] ) / mu_air_dry )

            # number density
            dict_atmprof['ndensity'] = dict_atmprof['plm'] / ( cgs.RR * dict_atmprof['TempL'] ) # mol

            # density
            layer_mu_atm = 1./ ( ( 1. - dict_atmprof['q'] ) / mu_air_dry + dict_atmprof['q'] / constants.MU_H2O )
            dict_atmprof['rho'] = dict_atmprof['ndensity'] * layer_mu_atm

            # refractivity
            refrac_layers  = molecules['H2O']['refractivity'] * ( dict_atmprof['ndensity'] * dict_atmprof['xH2O'] ) / cgs.amagat 
            if l_O3 :
                refrac_layers  = molecules['O3']['refractivity']  * ( dict_atmprof['ndensity'] * dict_atmprof['xO3'] ) / cgs.amagat 

            sum_layers     = deepcopy( dict_atmprof['xH2O'] )
            refrac_otherwise = 0.
            for molename in dict_NonCondensableGas :
                if dict_NonCondensableGas[molename]=='otherwise' :
                    refrac_otherwise = molecules[molename]['refractivity']
                else :
                    refrac_layers += molecules[molename]['refractivity'] * ( dict_atmprof['ndensity'] * ( 1. - dict_atmprof['xH2O'] ) * dict_NonCondensableGas[molename] ) / cgs.amagat
                    sum_layers    += ( 1. - dict_atmprof['xH2O'] ) * dict_NonCondensableGas[molename]
            refrac_layers += refrac_otherwise * ( dict_atmprof['ndensity'] * ( 1. - sum_layers ) ) / cgs.amagat
            dict_atmprof['refractivity'] = deepcopy( refrac_layers )

            # differential refractivity
            dndr_layers     = np.zeros_like( refrac_layers )
            dndr_layers[0]  = ( refrac_layers[1] - refrac_layers[0] ) / ( dict_atmprof['z'][1] - dict_atmprof['z'][0] )
            for zi in xrange( 1, len( refrac_layers )-1 ):
                dndr_layers[zi] = ( refrac_layers[zi+1] - refrac_layers[zi-1] ) / ( dict_atmprof['z'][zi+1] - dict_atmprof['z'][zi-1] )
            dndr_layers[-1] = ( refrac_layers[-1] - refrac_layers[-2] ) / ( dict_atmprof['z'][-1] - dict_atmprof['z'][-2] )
            dict_atmprof['dndr'] = -1.*deepcopy( dndr_layers ) # flip the sign

            # append
            list_dict_atmprof.append( dict_atmprof )

    return list_theta, list_dict_atmprof


#=============================================================================
def extrapolate_to_p_min( p_min, dict_atmprof, g_planet, param, mu_air_dry, n_toplayer=10 ):

    z_max      = dict_atmprof[   'z'][-1]
    q_max      = dict_atmprof[   'q'][-1]
    plm_max    = dict_atmprof[ 'plm'][-1]
    mu_atm_max = 1./ ( ( 1. - q_max ) / mu_air_dry + q_max / constants.MU_H2O )

    # additional layers for pressure
    scale_height_max = ( cgs.RR * dict_atmprof['TempL'][-1] ) / ( g_planet * mu_atm_max )
    plm_toplayers    = np.logspace( np.log10( plm_max ), np.log10( p_min ), n_toplayer+1  )[1:]

    dict_atmprof['plm'] = np.r_[ dict_atmprof['plm'], plm_toplayers ]

    # additional layers for other parameters
    for key in param :

        # altitude increases
        if key=='z' :
            z_toplayers       = z_max + scale_height_max * np.log( plm_max / plm_toplayers )
            dict_atmprof['z'] = np.r_[ dict_atmprof['z'], z_toplayers ]

        # otherwise parameters do not change
        else :
            add_toplayers     = np.tile( dict_atmprof[key][-1], n_toplayer )
            dict_atmprof[key] = np.r_[ dict_atmprof[key], add_toplayers ]

    return dict_atmprof



#=============================================================================
def truncate_to_p_min( p_min, dict_atmprof, param ):

    for key in param :
        dict_atmprof[key] = dict_atmprof[key][np.where( dict_atmprof['plm'] > p_min )]

    return dict_atmprof


