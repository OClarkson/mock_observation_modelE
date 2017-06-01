import numpy as np
import datetime
import os
import sys
import commands

from setup import *

import io_nc
import io_txt

from geometry import *
import plot

sec2hr = 1. / ( 60.*60. )

label_month = [ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' ]

#=======================================================================
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


    #-----------------------------------------------
    # Geometry set-up
    #-----------------------------------------------

    oblqty_deg, p_spin_sec, p_orbit_sec = io_txt.extract_param( RFILE, [ 'obliquity', 'siderealrotationperiod', 'siderealorbitalperiod' ], type='float' )

    if FULL_PHASE :
        PHASE0 = -180.
        TIME_END  = p_orbit_sec * sec2hr
        DT        = TIME_END / DIV_ORBIT

    oblqty, phase_eq, phase0, inc = init_geometry( oblqty_deg, PHASE_EQ, PHASE0, INC )
    omega_spin, omega_orbit       = init_omega( p_spin_sec * sec2hr, p_orbit_sec * sec2hr )

    datafile_example = commands.getoutput( "find " + DATAFILE_DIR[:-1] + " -name *"+DATAFILE_TAG ).split("\n")[0]
    nlat, nlon, array_lat, array_lon, array_data_dummy = io_nc.read_nc( datafile_example, mode="sw" )
    array_area  = init_area( nlat, nlon, array_lat, array_lon )
    lon_offset = get_lon_offset( array_lat, array_lon, SUBSTELLAR_LON_0, oblqty, phase_eq, phase0 )

    #-----------------------------------------------
    # Computing time variation...
    #-----------------------------------------------
    print 'Computing time variation...'
    #-----------------------------------------------
    
    data_integrated_sw = []
    data_integrated_lw = []

    list_time = []
    time      = 0.

    if MONTHLY :    
        month1_old = -1

    else :

        if ShortWave_LightCurve or ShortWave_Spectrum :
            nlat, nlon, array_lat, array_lon, array_data_sw = io_nc.read_nc( DATAFILE_DIR+'ANN'+DATAFILE_TAG, mode='sw' )

        if LongWave_LightCurve or LongWave_Spectrum :
            nlat, nlon, array_lat, array_lon, array_data_lw = io_nc.read_nc( DATAFILE_DIR+'ANN'+DATAFILE_TAG, mode='lw' )

    while time <= TIME_END :

        if MONTHLY :
            fraction_of_year = ( phase0 + omega_orbit*time - phase_eq ) / ( 2. * np.pi ) + ( ( 31+28+21 ) / 365.25 ) # Vernal equnox ~ March 21th
            month1 = int( np.floor( 12 * fraction_of_year - 0.5 ) )
            month2 = month1 + 1
            weight = 12 * fraction_of_year - 0.5 - month1
            month1 = month1 % 12
            month2 = month2 % 12

            if month1_old != month1 :
                print 'reading ' + label_month[month1]

                if ShortWave_LightCurve or ShortWave_Spectrum :
                    nlat, nlon, array_lat, array_lon, array_data_sw_1 = io_nc.read_nc( DATAFILE_DIR+label_month[month1]+DATAFILE_TAG, mode='sw' )
                    nlat, nlon, array_lat, array_lon, array_data_sw_2 = io_nc.read_nc( DATAFILE_DIR+label_month[month2]+DATAFILE_TAG, mode='sw' )
                    array_data_sw = array_data_sw_1 * ( 1. - weight ) + array_data_sw_2 * weight

                if LongWave_LightCurve or LongWave_Spectrum :
                    nlat, nlon, array_lat, array_lon, array_data_lw_1 = io_nc.read_nc( DATAFILE_DIR+label_month[month1]+DATAFILE_TAG, mode='lw' )
                    nlat, nlon, array_lat, array_lon, array_data_lw_2 = io_nc.read_nc( DATAFILE_DIR+label_month[month2]+DATAFILE_TAG, mode='lw' )
                    array_data_lw = array_data_lw_1 * ( 1. - weight ) + array_data_lw_2 * weight

                month1_old = month1

        #-----------------------------------------------
        # compute weight function
        #-----------------------------------------------
        array_cosTH0, array_cosTH1 = get_weight( omega_spin, omega_orbit, oblqty, phase_eq, phase0, inc, 
                                                 array_lat, array_lon + lon_offset, time )
        
        array_cosTH0[ np.where( array_cosTH0 < 0. )[0] ] = 0.
        array_cosTH1[ np.where( array_cosTH1 < 0. )[0] ] = 0.

        if ShortWave_LightCurve or ShortWave_Spectrum :
            array_weight = ( 1./np.pi ) * array_cosTH0 * array_cosTH1 * array_area
            data_integrated_sw.append( np.dot( array_weight, array_data_sw ) / np.sum( array_weight ) )

        if LongWave_LightCurve or LongWave_Spectrum :
            array_weight = ( 1./np.pi ) * array_cosTH1 * array_area
            data_integrated_lw.append( np.dot( array_weight, array_data_lw ) / np.pi )  # per area

        if FULL_PHASE :
            list_time.append( phase0 + omega_orbit*time )
        else :
            list_time.append( time )

        time = time + DT

    data_integrated_sw = np.array( data_integrated_sw )
    data_integrated_lw = np.array( data_integrated_lw )

    #-----------------------------------------------
    # plot
    #-----------------------------------------------

    if not DEBUG_ON :

        outfile_head = OUTFILE_DIR + OUTFILE_TAG

        if ShortWave_LightCurve :
            plot.plot_lc( list_time, data_integrated_sw, outfile_head, RFILE, SPECTRAL_FILE_SW, mode='sw', full_phase=FULL_PHASE )

        if ShortWave_Spectrum :
            plot.plot_sp( list_time, data_integrated_sw, outfile_head, RFILE, SPECTRAL_FILE_SW, mode='sw' )

        if LongWave_LightCurve :
            plot.plot_lc( list_time, data_integrated_lw, outfile_head, RFILE, SPECTRAL_FILE_LW, mode='lw', full_phase=FULL_PHASE )

        if LongWave_Spectrum :
            plot.plot_sp( list_time, data_integrated_lw, outfile_head, RFILE, SPECTRAL_FILE_LW, mode='lw' )


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

