#=============================================================================
# Module
#=============================================================================

LOOPMAX = 1000

import numpy as np
import multiprocessing
import functools
import datetime
import sys
import os
import matplotlib.pyplot as plt
from copy import deepcopy

from setup import *
from molecules import *
import transmission

import read_xs
import read_aijl
import read_ascii
import cgs
import resolution

import util_errors
import util_interp
import util_plot



#=============================================================================
# functions
#=============================================================================
def call_transmission( position_index, params ):

    grid_wn, list_theta, list_dict_atmprof, dict_griddata_logXSofWNTP, dict_NonCondensableGas, dict_geom = params
    theta        = list_theta[ position_index ]
    d_theta      = 2. * np.pi / ( len( list_theta ) )
    dict_atmprof = list_dict_atmprof[ position_index ]

    print 'Working on theta =', theta/np.pi*180.

    #------------------------------------------------
    # check range of lookup table
    #------------------------------------------------
    if MOLABS_ON :
        p_lookuptable_min = np.exp( np.min( dict_griddata_logXSofWNTP['coords'][:,1] ) )
        p_lookuptable_max = np.exp( np.max( dict_griddata_logXSofWNTP['coords'][:,1] ) )
        t_lookuptable_min = np.min( dict_griddata_logXSofWNTP['coords'][:,0] )
        if (dict_atmprof['plm'][0] > p_lookuptable_max ):
            util_errors.exit_msg( 'P'+str(position_index)+': maximum pressure of atmosphere is larger than maximum pressure in lookuptable.')
        if ( dict_atmprof['plm'][-1] < p_lookuptable_min ):
            util_errors.warning_longmsg( [ 'P'+str(position)+': minimum pressure of atmosphere is smaller than minimum pressure in lookuptable.' ,
                                           '    Cross section above '+str(p_lookuptable_min/cgs.mbar_to_barye)+' will be ignored.' ] )

    #------------------------------------------------
    # compute spectra
    #------------------------------------------------
    matrixW_Ftransmit = transmission.raytrace_opacity( grid_wn, theta, d_theta, dict_atmprof, dict_griddata_logXSofWNTP, DICT_NonCondensableGas, Z_TOP, Z_NUM, B_NUM, D_L, dict_geom, flag_molabs=MOLABS_ON, flag_rayleigh=RAYLEIGH_ON, flag_cld=CLD_ON, cld_d_eff_ice=CLD_D_EFF_ICE, cld_d_eff_liquid=CLD_D_EFF_LIQUID )

    return matrixW_Ftransmit


#=============================================================================
def call_multicore( list_index, params ):

    p = multiprocessing.Pool( CORE_NUM )
    result = p.map( functools.partial( call_transmission, params=params ), list_index )
    return result



#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":

    #------------------------------------------------
    # initialization 
    #------------------------------------------------

    print ''
    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")

    if not DEBUG_ON :

        # Create directory
        out_dir = OUTFILE_DIR + OUTFILE_TAG + "/"
        if os.path.exists( out_dir ):
            print out_dir+" already exists. Overwrite? [y/n] ...", 
            answer = raw_input()
            if answer=='n'  :
                sys.exit()
            elif not answer=='y' :
                util_errors.exit_msg('Unknown answer')
        else :
            os.mkdir( out_dir )
            print "Created directory:", out_dir

        # Save THIS file and the setup file for reproducibility
        # ( This idea comes from Jacob )
        thisfile  = os.path.basename(__file__)
        setupfile = "setup.py"
        os.system( "cp " + thisfile  + ' ' + out_dir + thisfile  )
        os.system( "cp " + setupfile + ' ' + out_dir + setupfile )
        print "Saved :", thisfile, " &", setupfile

        # Save start time
        filename_log = OUTFILE_DIR + OUTFILE_TAG + "/time.log"
        f_log = open( filename_log, 'w')
        f_log.write( "start time: " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n" )


    #---------------------------------------------------
    # GEOMETRICAL PARAMETERS
    #---------------------------------------------------

    dict_geom = { 'r_planet'        : R_PLANET, 
                  'r_star'          : R_STAR, 
                  'distance_to_star': DISTANCE_TO_STAR, 
                  'impact_parameter': IMPACT_PARAMETER, 
                  'orbital_phase'   : ORBITAL_PHASE }


    #------------------------------------------------
    # input atmospheric profiles
    #------------------------------------------------

    if '.nc' in DATAFILE :
        list_theta, list_dict_atmprof = read_aijl.extract_limbprof( DATAFILE, DICT_NonCondensableGas, Z_TOP, G_PLANET )
    elif '.txt' in DATAFILE : 
        list_theta, list_dict_atmprof = read_ascii.extract_prof( DATAFILE, DICT_NonCondensableGas, Z_TOP, G_PLANET )
    else :
        util_errors.exit_msg('Unknown format for DATAFILE. ')

    if ATMAVE_ON :

        if not ( len( list_theta ) % ATMAVE_NUM == 0 ) :
            util_errors.exit_msg('Irrelevant ATMAVE_NUM. ')

        list_ave_dict_atmprof = []
        list_ave_theta = []
        jj = 0
        for jj in xrange( len( list_dict_atmprof ) ):

            if ( jj % ATMAVE_NUM == 0 ) :
                list_ave_theta.append( deepcopy( list_theta[jj] ) )
                list_ave_dict_atmprof.append( deepcopy( list_dict_atmprof[jj] ) )
            else :
                list_ave_theta[-1] += list_theta[jj]
                for key in list_dict_atmprof[0] :
                    list_ave_dict_atmprof[-1][key] += list_dict_atmprof[jj][key]
            if ( jj % ATMAVE_NUM == ATMAVE_NUM - 1 ) :
                list_ave_theta[-1] = list_ave_theta[-1] / ( 1. * ATMAVE_NUM )
                for key in list_dict_atmprof[0] :
                    list_ave_dict_atmprof[-1][key] = list_ave_dict_atmprof[-1][key] / ( 1. * ATMAVE_NUM )
                
        list_dict_atmprof = deepcopy( list_ave_dict_atmprof )
        list_theta        = deepcopy( list_ave_theta )


    if not REFRACTION_ON :

        for ii in xrange( len( list_dict_atmprof ) ):
            list_dict_atmprof[ii]['dndr'] = np.zeros_like( list_dict_atmprof[0]['dndr'] )

    #------------------------------------------------
    # Input cross section tables for gases.
    # Set wavenumber grid.
    #------------------------------------------------

    if MOLABS_ON :

        # list of molecules
        list_mol = [ 'H2O', 'O3' ] + DICT_NonCondensableGas.keys()
#        list_mol = [ 'H2O' ] + DICT_NonCondensableGas.keys()

        # Adjust wavenumber range
        # if the wavenumber range of the lookup table is narrower, 
        # the wavenumber range of the lookup table is adopted
        tmp = np.load(  XSFILE_TAG + list_mol[0] + ".npz" )
        wn_min, wn_max, wn_num = tmp['WN'][0], tmp['WN'][-1], len(tmp['WN'])

        # get gridded data of cross section
        grid_wn, dict_griddata_logXSofWNTP = read_xs.griddata( list_mol, XSFILE_TAG, ( wn_min, wn_max ), cnt_h2o_on=CNTNM_H2O_ON )

    else :

        grid_wn = np.linspace( WN_MIN, WN_MAX, WN_NUM )
        dict_griddata_logXSofWNTP = {}
        


    #---------------------------------------------------
    # MAIN PART. 
    # calculate transmission at different locations
    #---------------------------------------------------

    params = ( grid_wn, list_theta, list_dict_atmprof, dict_griddata_logXSofWNTP, DICT_NonCondensableGas, dict_geom )
    list_index = range( len( list_dict_atmprof ) )
    
    if MULTICORE_ON :
        matrixW_Ftransmit = call_multicore( list_index, params )
        matrixW_Ftransmit_total = np.sum( np.array( matrixW_Ftransmit ), axis=0 )
        
    else :
        grid_sp_serial = np.zeros( len( list_index ) )
        matrixW_Ftransmit_total = np.zeros( len( grid_wn ) )
        for ii in list_index :
            matrixW_Ftransmit_total += call_transmission( ii, params )
            

    Fstar = np.pi
    planet_shadow = np.pi * ( dict_geom['r_planet'] + Z_TOP )**2 / ( dict_geom['r_star']**2 )
    matrixW_Ftransit = Fstar - planet_shadow + matrixW_Ftransmit_total

    matrixW_dF    = Fstar - matrixW_Ftransit 
    matrixW_Heff  = dict_geom['r_star'] * np.sqrt( matrixW_dF / np.pi ) - dict_geom['r_planet'] 
    matrixW_Heff  = matrixW_Heff  * 1e-5 # cm => km
    matrixW_dFppm = (matrixW_dF/Fstar)*1e6 # ppm


    if ( not DEBUG_ON ):
        data     = np.c_[ 1e4/grid_wn, matrixW_Heff, matrixW_dFppm  ]
        myheader = '# wavelength [um]\teffective altitude [km]\ttransit depth [ppm]'
        np.savetxt( OUTFILE_DIR + OUTFILE_TAG + '/transit_Heff_depth', data, header=myheader )

    #------------------------------------------------
    # Post-processing
    #------------------------------------------------

    # Lower resolution.
    if (LOW_RES_ON):
        grid2_wn, matrixW2_Heff  = resolution.lower_resolution( grid_wn, matrixW_Heff, RESOLUTION )
        grid2_wn, matrixW2_dFppm = resolution.lower_resolution( grid_wn, matrixW_dFppm, RESOLUTION )
        data2                    = np.c_[ 1e4/grid2_wn, matrixW2_Heff, matrixW2_dFppm ]
        if ( not DEBUG_ON ):
            np.savetxt( OUTFILE_DIR+OUTFILE_TAG+"/transit_Heff_depth_R"+str(RESOLUTION), data2, header=myheader )

    #------------------------------------------------
    # End.
    #------------------------------------------------
    # Save end time
    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")
    print ''

    if not DEBUG_ON :

        f_log.write( "end time: " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n" )
        f_log.close()

