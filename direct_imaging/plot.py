import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import errors
import sys
import io_txt

params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)


#=======================================================================
def get_band( spectral_file, mode, band_num ):

    band = io_txt.extract_block( spectral_file, 'Band        Lower limit         Upper limit', '*END' )
    if ( len( band ) != band_num ) :
        errors.exit_msg( 'Band number in spectral file does not match. ' )
    return band[:,1:] * 1e6 # m => um


#=======================================================================
def plot_sp( data_time, data_integrated, outfile_head, rfile, spectral_file, mode ):

    #-----------------------------------------------
    print 'Plotting average spectra'
    #-----------------------------------------------

    data_sp     = np.average( data_integrated, axis=0 )
    band_num    = len( data_sp )

    band        = get_band( spectral_file, mode, band_num ) 
    band_center = ( band[:,0] + band[:,1] ) * 0.5
    band_width  = ( band[:,1] - band[:,0] )

    fig, ax = plt.subplots(1,1,figsize=(5,3))
    ax.set_xlabel( 'wavelength [ $\mu $m ]' )

    if mode.lower()=='sw' :
        ax.set_xlim( [ 0.2, 2.0 ] )
        ax.set_ylabel( r'apparent albedo' )

    elif mode.lower()=='lw' :
        data_sp     = data_sp / band_width # per micron
        ax.set_xlim( [ 3., 20.0 ] )
#        ax.set_xscale( 'log' )
        ax.set_ylabel( r'radiation [W/m$^2$/$\mu$m/sr]' )

    ax.plot( band_center, data_sp, color='k' )

    plt.savefig( outfile_head+'/sp_'+mode+'.pdf', bbox_inches='tight' ) 

    data = np.c_[ band_center, data_sp ]
    np.savetxt( outfile_head+'/sp_'+mode+'.txt', data, header='wl\tsp' )


#=======================================================================
def plot_lc( data_time, data_integrated, outfile_head, rfile, spectral_file, mode, full_phase=True ):

    #-----------------------------------------------
    print 'Plotting light curves'
    #-----------------------------------------------

    band_num    = len( data_integrated.T )
    band        = get_band( spectral_file, mode, band_num ) 

    band        = get_band( spectral_file, mode, band_num ) 
    band_center = ( band[:,0] + band[:,1] ) * 0.5
    band_width  = ( band[:,1] - band[:,0] )

    fig, ax = plt.subplots(1,1,figsize=(5,3))
    if mode=='SW' or mode=='sw' :
        ax.set_ylabel( r'apparent albedo' )
        ax.set_ylim( [ 0., 1.0 ] )
    elif mode=='LW' or mode=='lw' :
        data_integrated = data_integrated / band_width # per micron
        ax.set_ylabel( r'radiation [W/m$^2$/$\mu$m/sr]' )
#        ax.set_ylabel( r'radiation' )

    if full_phase :
        ax.set_xlabel( 'orbital phase [rad]' )
        ax.set_xlim([-np.pi, np.pi])
    else :
        ax.set_xlabel( 'time [hr]' )

    for jj in xrange( len( data_integrated.T ) ):
        ax.plot( data_time, data_integrated.T[jj], label=str(band_center[jj])+r'$\mu $m' )

#    handles, labels = ax.get_legend_handles_labels()
#    ax.legend(handles, labels)

    plt.savefig( outfile_head+'/lc_'+mode+'.pdf', bbox_inches='tight' ) 

    data_time = np.array( data_time )
    data = np.c_[ data_time, data_integrated ]
    np.savetxt( outfile_head+'/lc_'+mode+'.txt', data, header='time\tsp' )
