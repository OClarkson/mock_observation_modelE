import numpy as np
import datetime
import os
import sys
import subprocess
from copy import deepcopy

from setup import *

import io_txt

import geometry
import plot

sec2hr = 1. / ( 60.*60. )

label_month = [ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' ]

if l_Socrates :
    from io_nc import read_nc_socrates as read_nc
else :
    from io_nc import read_nc_giss as read_nc

#=======================================================================
if __name__ == "__main__":

    #------------------------------------------------
    # initialization 
    #------------------------------------------------

    print('')
    now = datetime.datetime.now()
    print(now.strftime("%Y-%m-%d %H:%M:%S"))

    if not l_Debug :

        # Create directory
        out_dir = s_outFile_Dir + s_outFile_Tag + "/"
        if os.path.exists( out_dir ):
            print(out_dir+" already exists. Overwrite? [y/n] ...", end=' ') 
            answer = input()

            if answer=='n'  :
                sys.exit()

            elif not answer=='y' :
                util_errors.exit_msg('Unknown answer')
        else :
            os.mkdir( out_dir )
            print("Created directory:", out_dir)

        # Save THIS file and the setup file for reproducibility
        # ( This idea comes from Jacob )
        thisfile  = os.path.basename(__file__)
        setupfile = "setup.py"
        os.system( "cp " + thisfile  + ' ' + out_dir + thisfile  )
        os.system( "cp " + setupfile + ' ' + out_dir + setupfile )
        print("Saved :", thisfile, " &", setupfile)

        # Save start time
        filename_log = out_dir + "time.log"
        f_log = open( filename_log, 'w')
        f_log.write( "start time: " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n" )


    #-----------------------------------------------
    # Geometry set-up
    #-----------------------------------------------

    oblqty_deg, p_spin_sec, p_orbit_sec = io_txt.extract_param( s_RFile, [ 'obliquity', 'siderealrotationperiod', 'siderealorbitalperiod' ], type='float' )


    if l_FullOrbit :
        f_PhaseAngle_Initial_deg = -180. 
        f_TimeLimit_hr           = p_orbit_sec * sec2hr
        f_TimeInterval_hr        = f_TimeLimit_hr / i_DivideOrbit


    oblqty, phase_eq, phase0, inc = geometry.init_geometry( oblqty_deg, f_PhaseAngle_Equinox_deg, f_PhaseAngle_Initial_deg, f_InclinationAngle_deg )
    omega_spin, omega_orbit       = geometry.init_omega( p_spin_sec * sec2hr, p_orbit_sec * sec2hr )

    datafile_example = subprocess.getoutput( "find " + s_aijFile_Dir[:-1] + " -name *" + s_aijFile_Tag ).split("\n")[0]
    nlat, nlon, array_lat, array_lon, array_data_dummy = read_nc( datafile_example, mode="sw" )
    array_area = geometry.init_area( nlat, nlon, array_lat, array_lon )
    lon_offset = geometry.get_lon_offset( array_lat, array_lon, f_SubStellarLongitude_Initial_deg, oblqty, phase_eq, phase0 )

    #-----------------------------------------------
    # Computing time variation...
    #-----------------------------------------------
    print('Computing time variation...')
    #-----------------------------------------------
    
    data_integrated_sw = []
    data_integrated_lw = []

    list_time = []
    time      = 0.

    if l_Monthly :    
        month1_old = -1

    else :
        ANNaijfile = s_aijFile_Dir + s_aijFile_Tag
        if l_ShortWave_LightCurve or l_ShortWave_Spectrum :
            nlat, nlon, array_lat, array_lon, array_data_sw = read_nc( ANNaijfile, mode='sw' )

        if l_LongWave_LightCurve or l_LongWave_Spectrum :
            nlat, nlon, array_lat, array_lon, array_data_lw = read_nc( ANNaijfile, mode='lw' )

    while time <= f_TimeLimit_hr :

        if l_Monthly :

            fraction_of_year = ( phase0 + omega_orbit*time - phase_eq ) / ( 2. * np.pi ) + ( ( 31+28+21 ) / 365.25 ) # Vernal equnox ~ March 21th
            month1 = int( np.floor( 12 * fraction_of_year - 0.5 ) ) # assuming that the monthly outputs correspond to the middle of the month
            month2 = month1 + 1
            weight = 12 * fraction_of_year - 0.5 - month1
            month1 = month1 % 12
            month2 = month2 % 12

            if month1_old != month1 :

                if l_ShortWave_LightCurve or l_ShortWave_Spectrum :
                    if month1_old == -1 :
                        nlat, nlon, array_lat, array_lon, array_data_sw_1 = read_nc( s_aijFile_Dir+label_month[month1]+s_aijFile_Tag, mode='sw' )
                        print('Reading ', label_month[month1], 'data for SW')
                    else :
                        array_data_sw_1 = deepcopy( array_data_sw_2 )
                    print('Reading ', label_month[month2], 'data for SW')
                    nlat, nlon, array_lat, array_lon, array_data_sw_2 = read_nc( s_aijFile_Dir+label_month[month2]+s_aijFile_Tag, mode='sw' )

                if l_LongWave_LightCurve or l_LongWave_Spectrum :
                    if month1_old == -1 :
                        nlat, nlon, array_lat, array_lon, array_data_lw_1 = read_nc( s_aijFile_Dir+label_month[month1]+s_aijFile_Tag, mode='lw' )
                        print('Reading ', label_month[month1], 'data for LW')
                    else :
                        array_data_lw_1 = deepcopy( array_data_lw_2 )
                    print('Reading ', label_month[month2], 'data for LW')
                    nlat, nlon, array_lat, array_lon, array_data_lw_2 = read_nc( s_aijFile_Dir+label_month[month2]+s_aijFile_Tag, mode='lw' )

                month1_old = month1

            if l_ShortWave_LightCurve or l_ShortWave_Spectrum :
                array_data_sw = array_data_sw_1 * ( 1. - weight ) + array_data_sw_2 * weight
                array_data_sw[ np.where( array_data_sw_1 == 0. )] = array_data_sw_2[ np.where( array_data_sw_1 == 0. )]
                array_data_sw[ np.where( array_data_sw_2 == 0. )] = array_data_sw_1[ np.where( array_data_sw_2 == 0. )]

            if l_LongWave_LightCurve or l_LongWave_Spectrum :
                array_data_lw = array_data_lw_1 * ( 1. - weight ) + array_data_lw_2 * weight

        #-----------------------------------------------
        # compute weight function
        #-----------------------------------------------
        array_cosTH0, array_cosTH1 = geometry.get_weight( omega_spin, omega_orbit, oblqty, phase_eq, phase0, inc, 
                                                          array_lat, array_lon + lon_offset, time )
        
        array_cosTH0[ np.where( array_cosTH0 < 0. )[0] ] = 0.
        array_cosTH1[ np.where( array_cosTH1 < 0. )[0] ] = 0.

        if l_ShortWave_LightCurve or l_ShortWave_Spectrum :

            array_weight = ( 1./np.pi ) * array_cosTH0 * array_cosTH1 * array_area

            if np.sum( array_weight ) == 0. :
                data_integrated_sw.append( np.zeros( len( array_data_sw[0] ) ) )
            else :
                data_integrated_sw.append( np.dot( array_weight, array_data_sw ) / np.sum( array_weight ) )

        if l_LongWave_LightCurve or l_LongWave_Spectrum :
            array_weight = ( 1./np.pi ) * array_cosTH1 * array_area
            data_integrated_lw.append( np.dot( array_weight, array_data_lw ) / np.pi )  # per area

        if l_FullOrbit :
            list_time.append( phase0 + omega_orbit*time )
        else :
            list_time.append( time )

        time = time + f_TimeInterval_hr

    data_integrated_sw = np.array( data_integrated_sw )
    data_integrated_lw = np.array( data_integrated_lw )

    #-----------------------------------------------
    # plot
    #-----------------------------------------------

    if not l_Debug :

        if l_ShortWave_LightCurve :
            if l_Socrates :
                plot.plot_lc( list_time, data_integrated_sw, out_dir, s_RFile, spfile=s_SpectralFile_SW, mode='sw', full_phase=l_FullOrbit )
            else :
                plot.plot_lc( list_time, data_integrated_sw, out_dir, s_RFile, spfile=False, mode='sw', full_phase=l_FullOrbit )

        if l_ShortWave_Spectrum :
            if l_Socrates :
                plot.plot_sp( list_time, data_integrated_sw, out_dir, s_RFile, spfile=s_SpectralFile_SW, mode='sw' )
            else :
                util_errors.exit_msg('Spectrum plot is not supported for giss radiation.')

        if l_LongWave_LightCurve :
            if l_Socrates :
                plot.plot_lc( list_time, data_integrated_lw, out_dir, s_RFile, spfile=s_SpectralFile_LW, mode='lw', full_phase=l_FullOrbit )
            else :
                plot.plot_lc( list_time, data_integrated_lw, out_dir, s_RFile, spfile=False, mode='lw', full_phase=l_FullOrbit )

        if l_LongWave_Spectrum :
            if l_Socrates :
                plot.plot_sp( list_time, data_integrated_lw, out_dir, s_RFile, spfile=s_SpectralFile_LW, mode='lw' )
            else :
                util_errors.exit_msg('Spectrum plot is not supported for giss radiation.')

    #------------------------------------------------
    # End.
    #------------------------------------------------
    # Save end time
    now = datetime.datetime.now()
    print(now.strftime("%Y-%m-%d %H:%M:%S"))
    print('')

    if not l_Debug :

        f_log.write( "end time: " + now.strftime("%Y-%m-%d %H:%M:%S") + "\n" )
        f_log.close()

