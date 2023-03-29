from setup import s_fig_format
import util_errors
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import setup_plot 

#=======================================================================
def plot_sp( outfile_head, dict_geom, grid_wl, matrixW_Heff, matrixW_dFppm ):

    #-----------------------------------------------
    print('Plotting transmission spectrum...')
    #-----------------------------------------------
    fig, ax, linecolor, title, y_unit, y_range = setup_plot.sp( grid_wl, matrixW_Heff, matrixW_dFppm )

    ax.set_ylim( y_range )

    if y_unit[0] == 'Heff' :
        ax.set_ylabel('effective radius [km]')
        ax.plot( grid_wl, matrixW_Heff, color=linecolor )

    elif y_unit[0] == 'ppm' :
        ax.set_ylabel('transit depth [ppm]')
        ax.plot( grid_wl, matrixW_dFppm, color=linecolor )

    if len( y_unit ) > 1 :
        ax2 = ax.twinx()
        outfile = outfile_head + '/sp_' + y_unit[0] + '_' + y_unit[1] + '.' + s_fig_format
        if y_unit[0] == 'Heff' and y_unit[1] == 'ppm' :
            ax2.set_ylabel('transit depth [ppm]')
            ppm_min = ( dict_geom['r_planet'] + y_range[0]*1e5 )**2 / ( dict_geom['r_star']**2 ) * 1e6
            ppm_max = ( dict_geom['r_planet'] + y_range[1]*1e5 )**2 / ( dict_geom['r_star']**2 ) * 1e6
            ax2.set_ylim([ ppm_min, ppm_max ])

        elif y_unit[0] == 'ppm' and y_unit[1] == 'Heff' :
            ax2.set_ylabel('effective radius [km]')
            Heff_min = np.sqrt( dict_geom['r_star']**2 * y_range[0]*1e-6 ) - dict_geom['r_planet'] 
            Heff_max = np.sqrt( dict_geom['r_star']**2 * y_range[1]*1e-6 ) - dict_geom['r_planet'] 
            ax2.set_ylim([ Heff_min, Heff_max ])

        else :
            util_errors.exit_msg( 'Inappropriate unit for y-axis.' )

    else :
        outfile = outfile_head + '/sp_' + y_unit[0] + '.' + s_fig_format

    if title :
        plt.title( title )
        
    plt.savefig( outfile, bbox_inches='tight' ) 

