import numpy as np
import datetime
import os
import io_nc
import io_txt
from geometry import *
from setup import *
import sys 
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

    nlat, nlon, array_lat, array_lon, array_data = io_nc.read_nc( DATAFILE_DIR+'JAN'+DATAFILE_TAG, mode=MODE )
    oblqty_deg, p_spin_sec, p_orbit_sec          = io_txt.extract_param( RFILE, [ 'OBLIQUITY', 'siderealRotationPeriod', 'siderealOrbitalPeriod' ], type='float' )
    oblqty, phase_eq, phase0, inc                = init_geometry( oblqty_deg, PHASE_EQ, PHASE0, INC )
    omega_spin, omega_orbit                      = init_omega( p_spin_sec * sec2hr, p_orbit_sec * sec2hr )
    array_area                                   = init_area( nlat, nlon, array_lat, array_lon )

    #-----------------------------------------------
    # Computing time variation...
    #-----------------------------------------------
    print 'Computing time variation...'
    #-----------------------------------------------

    data_integrated = []
    list_time       = []
    time = 0.
    
    month1_old = -1
    while time <= TIME_END :

        fraction_of_year = ( phase0 + omega_orbit*time - phase_eq ) / ( 2. * np.pi ) + ( ( 31+28+21 ) / 365.25 ) # Vernal equnox ~ March 21th
        month1 = int( np.floor( 12 * fraction_of_year - 0.5 ) )
        month2 = month1 + 1
        weight = 12 * fraction_of_year - 0.5 - month1
        month1 = month1 % 12
        month2 = month2 % 12
        if month1_old != month1 :
            print 'reading ' + label_month[month1]
            nlat, nlon, array_lat, array_lon, array_data1 = io_nc.read_nc( DATAFILE_DIR+label_month[month1]+DATAFILE_TAG, mode=MODE )
            nlat, nlon, array_lat, array_lon, array_data2 = io_nc.read_nc( DATAFILE_DIR+label_month[month2]+DATAFILE_TAG, mode=MODE )
            month1_old = month1

        array_data = array_data1 * ( 1. - weight ) + array_data2 * weight

        #-----------------------------------------------
        # compute weight function
        #-----------------------------------------------
        array_cosTH0, array_cosTH1 = get_weight( omega_spin, omega_orbit, oblqty, phase_eq, phase0, inc, 
                                                 array_lat, array_lon, time )
        
        array_cosTH0[ np.where( array_cosTH0 < 0. )[0] ] = 0.
        array_cosTH1[ np.where( array_cosTH1 < 0. )[0] ] = 0.

        if ( MODE == 'SW' or MODE == 'sw' ) :
            array_weight = ( 1./np.pi ) * array_cosTH0 * array_cosTH1 * array_area
            data_integrated.append( np.dot( array_weight, array_data ) / np.sum( array_weight ) )

        elif ( MODE == 'LW' or MODE == 'lw' ) :
            array_weight = ( 1./np.pi ) * array_cosTH1 * array_area
            data_integrated.append( np.dot( array_weight, array_data ) / np.pi )  # per area

        list_time.append( time )
        time = time + DT

    data_integrated = np.array( data_integrated )

    #-----------------------------------------------
    # plot
    #-----------------------------------------------

    if not DEBUG_ON :
        outfile_head = OUTFILE_DIR + OUTFILE_TAG
        if PLOT_LC :
            plot.plot_lc( list_time, data_integrated, outfile_head, RFILE, SPECTRAL_FILE, MODE )
        if PLOT_SP :
            plot.plot_sp( list_time, data_integrated, outfile_head, RFILE, SPECTRAL_FILE, MODE )


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

