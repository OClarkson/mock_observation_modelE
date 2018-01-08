import numpy as np
import util_errors


#=============================================================================
def gaussian( wl, wl0, resolution ):

    FWHM   = wl0 / resolution
    sigma  = FWHM / ( 2. * np.sqrt( 2. * np.log( 2 ) ) )
    result = 1. / np.sqrt( 2. * np.pi * sigma**2 ) * np.exp( - ( wl - wl0 )**2 / ( 2. * sigma**2 ) )
    return result


#=============================================================================
def lower_resolution( grid_wl, matrixW_sp, resolution ):

    # integral[ exp( - ( wl_i - wl_center) / sigma_center ) * sp ( wl_i ) * d wl_i ]

    mesh_wl, mesh_wl0 = np.meshgrid( grid_wl, grid_wl, indexing='ij' )
    matrixW_convfunc  = gaussian( mesh_wl, mesh_wl0, resolution )

    matrixW_central = np.r_[ grid_wl[0], 0.5 * ( grid_wl[:-1] + grid_wl[1:] ), grid_wl[-1] ]
    matrixW_dwl     = np.diff( matrixW_central )
    matrixW_sp_new  = np.dot( matrixW_convfunc, matrixW_sp * matrixW_dwl ) 

    # FWHM > 0.5 * delta_wl ?
    list_index = np.where( grid_wl[1:-1]/resolution < 0.5 * matrixW_dwl[1:-1] )[0]
    if len( list_index ) > 0 :
        util_errors.warning_longmsg( ['Due to the coarser grids of look-up tables than specified resolution, ',
                                      'error in low-resolution spectrum is large beyond ~' + str( int( grid_wl[list_index[0]] ) ) + 'um.', ] )

    return matrixW_sp_new


