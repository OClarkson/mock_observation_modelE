#------------------------------------------------
# ShortWave Spectrum
#------------------------------------------------
def sp( ):

    # figure size
    fig, ax = plt.subplots(1,1,figsize=(5,3))


    # color 
    linecolor = 'black'

    # title
    title = False

    # x axis
    ax.set_xlim( [ 1.0, 10. ] ) # 
    ax.set_xlabel( r'wavelength [$\mu $m]' )

    # y axis
    y_unit = [ 'Heff', 'ppm' ] 
    y_range = [ 0., 60. ]

    return fig, ax, linecolor, title, y_unit, y_range

