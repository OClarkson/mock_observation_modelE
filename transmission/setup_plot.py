import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------
# ShortWave Spectrum
#------------------------------------------------
def sp( grid_wl, matrixW_Heff, matrixW_dFppm ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))

    print(matrixW_Heff)
    # color 
    linecolor = 'black'

    # title
    title = False

    # x axis
    wl_max = np.minimum( grid_wl[-1], 10. )
    ax.set_xlim( [ grid_wl[0], wl_max ] ) # 
    ax.set_xlabel( r'wavelength [$\mu $m]' )

    # Uncomment the line below for logscale x-axis
    # ax.set_xscale('log')
    
    # y axis
    y_unit = [ 'Heff', 'ppm' ] 
    if y_unit[0] == 'Heff' :
        Heff_max = np.max( matrixW_Heff[np.where( grid_wl < wl_max )] ) + 10. 
        y_range = [ 0., Heff_max ]
    elif y_unit[0] == 'ppm' :
        dFppm_min = np.min( matrix_dFppm )
        dFppm_max = np.max( matrix_dFppm )
        dFppm_range = dFppm_max - dFppm_min
        y_range = [ dFppm_min - dFppm_range / 3.,  dFppm_max + dFppm_range / 3. ]

    return fig, ax, linecolor, title, y_unit, y_range

