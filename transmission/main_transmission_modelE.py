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
import plot

import util_errors
import util_interp

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
    if l_molecular_absorption :
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
    matrixW_Ftransmit = transmission.raytrace_opacity( grid_wn, theta, d_theta, dict_atmprof, dict_griddata_logXSofWNTP, DICT_NonCondensableGas, 
                                                       f_Z_top, i_Z_num, i_B_num, f_dL, dict_geom, 
                                                       flag_molabs=l_molecular_absorption, 
                                                       flag_rayleigh=l_Rayleigh, 
                                                       flag_cld=l_cloud )
#, 
#                                                       cld_d_eff_liquid=f_cloud_Deff_liquid, 
#                                                       cld_d_eff_ice=f_cloud_Deff_ice )

    return matrixW_Ftransmit


#=============================================================================
def call_multicore( list_index, params ):

    p = multiprocessing.Pool( i_core_num )
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

    if not l_Debug :

        # Create directory
        out_dir = s_outFile_Dir + s_outFile_Tag + "/"
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
        filename_log = s_outFile_Dir + s_outFile_Tag + "/time.log"
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

    if '.nc' in s_atmFile :
        list_theta, list_dict_atmprof = read_aijl.extract_limbprof( s_atmFile, DICT_NonCondensableGas, f_Z_top, G_PLANET )
    elif '.txt' in s_atmFile : 
        list_theta, list_dict_atmprof = read_ascii.extract_prof( s_atmFile, DICT_NonCondensableGas, f_Z_top, G_PLANET )
    else :
        util_errors.exit_msg('Unknown format for DATAFILE. ')

    if l_atm_average :

        if not ( len( list_theta ) % i_atmave_num == 0 ) :
            util_errors.exit_msg('Irrelevant i_atmave_num. ')

        list_ave_dict_atmprof = []
        list_ave_theta = []
        jj = 0
        for jj in xrange( len( list_dict_atmprof ) ):

            if ( jj % i_atmave_num == 0 ) :
                list_ave_theta.append( deepcopy( list_theta[jj] ) )
                list_ave_dict_atmprof.append( deepcopy( list_dict_atmprof[jj] ) )
            else :
                list_ave_theta[-1] += list_theta[jj]
                for key in list_dict_atmprof[0] :
                    list_ave_dict_atmprof[-1][key] += list_dict_atmprof[jj][key]
            if ( jj % i_atmave_num == i_atmave_num - 1 ) :
                list_ave_theta[-1] = list_ave_theta[-1] / ( 1. * i_atmave_num )
                for key in list_dict_atmprof[0] :
                    list_ave_dict_atmprof[-1][key] = list_ave_dict_atmprof[-1][key] / ( 1. * i_atmave_num )
                
        list_dict_atmprof = deepcopy( list_ave_dict_atmprof )
        list_theta        = deepcopy( list_ave_theta )


    if not l_refraction :

        for ii in xrange( len( list_dict_atmprof ) ):
            list_dict_atmprof[ii]['dndr'] = np.zeros_like( list_dict_atmprof[0]['dndr'] )

    #------------------------------------------------
    # Input cross section tables for gases.
    # Set wavenumber grid.
    #------------------------------------------------

    if l_molecular_absorption :

        # list of molecules
        if l_O3 :
            list_mol = [ 'H2O', 'O3' ] + DICT_NonCondensableGas.keys()
        else :
            list_mol = [ 'H2O' ] + DICT_NonCondensableGas.keys()

        # Adjust wavenumber range
        # if the wavenumber range of the lookup table is narrower, 
        # the wavenumber range of the lookup table is adopted
        for mol in list_mol :
            xsfile_tmp = s_xsFile_Tag + mol + ".npz" 
            if os.path.exists( xsfile_tmp ) : 
                break
        tmp = np.load( xsfile_tmp )

        grid_wn_org = tmp['WN']
        grid_T      = tmp['T']
        grid_P      = tmp['P']*cgs.mbar_to_barye

        indx_min = read_xs.nearestindex_WN( grid_wn_org, f_Wavenumber_min )
        indx_max = read_xs.nearestindex_WN( grid_wn_org, f_Wavenumber_max  )
        if indx_min < indx_max :
            grid_wn  = grid_wn_org[indx_min:indx_max+1]
        else :
            util_errors.exit_msg( 'Specified wavenumber range is not valid---check the consistensy with the cross section table. ')

        # get gridded data of cross section from HITRAN
        dict_griddata_logXSofWNTP = read_xs.griddata_line( list_mol, s_xsFile_Tag, grid_wn, grid_T, grid_P, cnt_h2o_on=l_H2O_continuum )

        if l_O3 and ( grid_wn[-1] > 1.e4 ) : # shortward of wavelength=1um
            dict_griddata_logXSofWNTP = read_xs.griddata_add_UV( 'O3', 'data/SerdyuchenkoGorshelevVersionJuly2013.dat', grid_wn,  grid_T, grid_P, dict_griddata_logXSofWNTP )

    else :

        grid_wn = np.linspace( f_Wavenumber_min, f_Wavenumber_max, i_Wavenumber_num )
        dict_griddata_logXSofWNTP = {}
        
    #---------------------------------------------------
    # MAIN PART. 
    # calculate transmission at different locations
    #---------------------------------------------------

    params = ( grid_wn, list_theta, list_dict_atmprof, dict_griddata_logXSofWNTP, DICT_NonCondensableGas, dict_geom )
    list_index = range( len( list_dict_atmprof ) )
    
    if l_multicore and ( i_core_num <= len( list_index ) ) :
        matrixW_Ftransmit = call_multicore( list_index, params )
        matrixW_Ftransmit_total = np.sum( np.array( matrixW_Ftransmit ), axis=0 )
        
    else :
        grid_sp_serial = np.zeros( len( list_index ) )
        matrixW_Ftransmit_total = np.zeros( len( grid_wn ) )
        for ii in list_index :
            matrixW_Ftransmit_total += call_transmission( ii, params )
            

    Fstar = np.pi
    planet_shadow = np.pi * ( dict_geom['r_planet'] + f_Z_top )**2 / ( dict_geom['r_star']**2 )
    matrixW_Ftransit = Fstar - planet_shadow + matrixW_Ftransmit_total

    matrixW_dF    = Fstar - matrixW_Ftransit 
    matrixW_Heff  = dict_geom['r_star'] * np.sqrt( matrixW_dF / np.pi ) - dict_geom['r_planet'] 
    matrixW_Heff  = matrixW_Heff  * 1e-5 # cm => km
    matrixW_dFppm = (matrixW_dF/Fstar)*1e6 # ppm

    # reverse the order
    grid_wl    = 1e4 / grid_wn[::-1]
    matrixW_Heff  = matrixW_Heff[::-1]
    matrixW_dFppm = matrixW_dFppm[::-1]

    #------------------------------------------------
    # Output
    #------------------------------------------------

    if ( not l_Debug ):

        data     = np.c_[ grid_wl, matrixW_Heff, matrixW_dFppm  ]
        myheader = '# wavelength [um]\teffective altitude [km]\ttransit depth [ppm]'
        np.savetxt( s_outFile_Dir + s_outFile_Tag + '/transit_Heff_depth', data, header=myheader )

        #------------------------------------------------
        # Post-processing
        #------------------------------------------------

        # Lower resolution.
        if l_lower_resolution :

            matrixW_wl, matrixW_Heff  = resolution.lower_resolution( grid_wl,  matrixW_Heff, f_resolution )
            matrixW_wl, matrixW_dFppm = resolution.lower_resolution( grid_wl, matrixW_dFppm, f_resolution )

            data  = np.c_[ matrixW_wl, matrixW_Heff, matrixW_dFppm ]
            np.savetxt( s_outFile_Dir + s_outFile_Tag + "/transit_Heff_depth_R"+str( int( f_resolution ) ), data, header=myheader )

            grid_wl = matrixW_wl

        # Make a plot
        if l_Plot :
            plot.plot_sp( out_dir, dict_geom, grid_wl, matrixW_Heff, matrixW_dFppm )

    #------------------------------------------------
    # End.
    #------------------------------------------------
    # Save end time
    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M:%S")
    print ''

    if not l_Debug :

        f_log.write( "end time: " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n" )
        f_log.close()

