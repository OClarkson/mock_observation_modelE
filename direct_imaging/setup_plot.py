import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------
# ShortWave Spectrum
#------------------------------------------------
def sp_sw( ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))

    # range
    ax.set_xlim( [ 0.3, 2.0 ] ) # 
    ax.set_ylim( [ 0.0, 1.0 ] )

    # label
    ax.set_xlabel( r'wavelength [$\mu $m]' )
    ax.set_ylabel( r'apparent albedo' )

    return fig, ax


#------------------------------------------------
# LongWave Spectrum
#------------------------------------------------
def sp_lw( ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))

    # range
    ax.set_xlim( [ 3., 20.0 ] )
    ax.set_ylim( [ 0., 10.0 ] )

    # label
    ax.set_xlabel( r'wavelength [$\mu $m]' )
    ax.set_ylabel( r'radiation [W/m$^2$/$\mu$m/sr]' )

    return fig, ax


#------------------------------------------------
# ShortWave Lightcurve
#------------------------------------------------
def lc_sw( data_num ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))

    # range
    ax.set_xlim( [ 0.3, 2.0 ] )
#    ax.set_ylim( [ 0.0, 1.0 ] )

    # labels
    ax.set_xlabel( r'time [hour]' )
    ax.set_ylabel( r'apparent albedo' )

    # colors
    if data_num > 1 :
        cmap   = plt.get_cmap( "rainbow" )
        colors = cmap( np.linspace( 0., 1., data_num ) )
    else :
        colors = 'k'

    return fig, ax, colors

#------------------------------------------------
def lc_sw_legend( ax ):

    handles, labels = ax.get_legend_handles_labels()
    ax.legend( handles, labels, 
               loc  = 3, # location
               ncol = 1, # number of columns
               prop = {'size'   :  8 ,
                       'family' : 'sans-serif', 
                       'style'  : 'normal'} # font properties
               )
    return ax


#------------------------------------------------
# LongWave Lightcurve
#------------------------------------------------
def lc_lw( data_num ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))

    # range
    ax.set_xlim( [ 3., 20. ] )
#    ax.set_ylim( [ 0., 6. ] )

    # label
    ax.set_xlabel( r'time [hour]' )
    ax.set_ylabel( r'radiation [W/m$^2$/$\mu$m/sr]' )

    # colors
    if data_num > 1 :
        cmap   = plt.get_cmap( "rainbow" )
        colors = cmap( np.linspace( 0., 1., data_num ) )
    else :
        colors = 'k'

    return fig, ax, colors

#------------------------------------------------
def lc_lw_legend( ax ):

    handles, labels = ax.get_legend_handles_labels()
    ax.legend( handles, labels, 
               loc  = 3, # location
               ncol = 1, # number of columns
               prop = {'size'   :  8 ,
                       'family' : 'sans-serif', 
                       'style'  : 'normal'} # font properties
               )
    return ax



