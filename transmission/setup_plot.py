import matplotlib.pyplot as plt
import numpy as np

#------------------------------------------------
# ShortWave Spectrum
#------------------------------------------------
def sp( grid_wl, matrixW_Heff, matrixW_dFppm ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))


    # color 
    linecolor = 'black'

    # title
    title = False

    # x axis
    ax.set_xlim( [ grid_wl[0], grid_wl[-1] ] ) # 
    ax.set_xlabel( r'wavelength [$\mu $m]' )

    # y axis
    y_unit = [ 'Heff', 'ppm' ] 
    if y_unit[0] == 'Heff' :
        y_range = [ 0., np.max( matrixW_Heff ) + 10. ]
    elif y_unit[0] == 'ppm' :
        dFppm_min = np.min( matrix_dFppm )
        dFppm_max = np.max( matrix_dFppm )
        dFppm_range = dFppm_max - dFppm_min
        y_range = [ dFppm_min - dFppm_range / 3.,  dFppm_max + dFppm_range / 3. ]

    return fig, ax, linecolor, title, y_unit, y_range

