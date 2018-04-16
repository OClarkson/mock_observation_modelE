import numpy as np
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import errors
import sys
import io_txt
import setup_plot

params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)


#=======================================================================
def get_band( spectral_file, mode, band_num ):

    band = io_txt.extract_block( spectral_file, "*BLOCK: TYPE =    1", '*END' )
    if ( len( band ) != band_num ) :
        errors.exit_msg( 'Band number in spectral file does not match. ' )
    return band[:,1:] * 1e6 # m => um



#=======================================================================
def get_labels( band, digit=0 ) :

    band       = np.array( band )
    band_shape = band.shape
    band_flat  = band.flatten()

    if digit == 0 :
        band_str_flat = map( str, map( int, map( round, band_flat ) ) )
    else :
        band_str_flat = map( str, map( lambda x: round(x, digit), band_flat ) )

    band_str_flat = np.array( band_str_flat )
    band_str      = band_str_flat.reshape( band_shape )

    labels = []
    for jj in xrange( len( band ) ):
        labels.append( band_str[jj][0] + '-' + band_str[jj][1] + r' $\mu $m' )

    return labels



#=======================================================================
def plot_sp( data_time, data_integrated, outfile_head, rfile, spfile=False, mode=False ):

    #-----------------------------------------------
    print 'Plotting ' + str(mode.upper()) +' average spectra' 
    #-----------------------------------------------

    data_sp     = np.nanmean( data_integrated, axis=0 )
    band_num    = len( data_sp )

    band        = get_band( spfile, mode, band_num ) 
    band_center = ( band[:,0] + band[:,1] ) * 0.5
    band_width  = ( band[:,1] - band[:,0] )

    if mode.lower()=='sw' :
        fig, ax = setup_plot.sp_sw( )

    elif mode.lower()=='lw' :
        fig, ax = setup_plot.sp_lw( )
        data_sp     = data_sp / band_width # per micron

    ax.plot( band_center, data_sp, color='k' )
    plt.savefig( outfile_head+'/sp_'+mode+'.pdf', bbox_inches='tight' ) 

    data = np.c_[ band_center, data_sp ]
    np.savetxt( outfile_head+'/sp_'+mode+'.txt', data, header='wl\tsp' )


#=======================================================================
def plot_lc( data_time, data_integrated, outfile_head, rfile, spfile=False, mode=False, full_phase=True ):

    #-----------------------------------------------
    print 'Plotting ' + str(mode.upper()) +' lightcurves' 
    #-----------------------------------------------

    band_num    = len( data_integrated.T )

    if band_num > 1 :
        band        = get_band( spfile, mode, band_num ) 
        band_center = ( band[:,0] + band[:,1] ) * 0.5
        band_width  = ( band[:,1] - band[:,0] )

    if mode=='SW' or mode=='sw' :
        fig, ax, colors = setup_plot.lc_sw( len( data_integrated.T ) )
        if band_num > 1 :
            labels  = get_labels( band, digit=2 )
    elif mode=='LW' or mode=='lw' :
        fig, ax, colors = setup_plot.lc_lw( len( data_integrated.T ) )
        if band_num > 1 :
            labels  = get_labels( band, digit=1 )
            data_integrated = data_integrated / band_width # per micron

    if full_phase :
        ax.set_xlabel( 'orbital phase [rad]' )
        ax.set_xlim([-np.pi, np.pi])

    for jj in xrange( len( data_integrated.T ) ):
        if mode=='SW' or mode=='sw' :
            if band_num > 1 :
                ax.plot( data_time, data_integrated.T[jj], label=labels[jj], c=colors[jj] )
            else :
                ax.plot( data_time, data_integrated.T[0], c=colors )

        elif mode=='LW' or mode=='lw' :
            if band_num > 1 :
                jj2 = len(data_integrated.T) - jj - 1
                ax.plot( data_time, data_integrated.T[jj2], label=labels[jj2], c=colors[jj] )
            else :
                ax.plot( data_time, data_integrated.T[0], c=colors )
    # legend
    if band_num > 1 :
        if mode=='SW' or mode=='sw' :
            ax = setup_plot.lc_sw_legend( ax )

        elif mode=='LW' or mode=='lw' :
            ax = setup_plot.lc_lw_legend( ax )

    plt.savefig( outfile_head+'/lc_'+mode+'.pdf', bbox_inches='tight' ) 

    data_time = np.array( data_time )
    data = np.c_[ data_time, data_integrated ]
    np.savetxt( outfile_head+'/lc_'+mode+'.txt', data, header='time\tsp' )
