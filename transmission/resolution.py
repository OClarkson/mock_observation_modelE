import numpy as np



#=============================================================================
def gaussian( wl, wl0, resolution ):

    FWHM   = wl0 / resolution
    sigma  = FWHM / ( 2. * np.sqrt( 2. * np.log( 2 ) ) )
    result = 1. / np.sqrt( 2. * np.pi * sigma**2 ) * np.exp( - ( wl - wl0 )**2 / ( 2. * sigma**2 ) )
    return result


#=============================================================================
def lower_resolution( grid_wl, matrixW_sp, resolution ):

    # integral[ exp( - ( wl_i - wl_center) / sigma_center ) * sp ( wl_i ) * d wl_i ]

    mesh_wl0, mesh_wl = np.meshgrid( grid_wl, grid_wl, indexing='ij' )
    matrixW_convfunc  = gaussian( mesh_wl, mesh_wl0, resolution )

    matrixW_central = np.r_[ grid_wl[0], 0.5 * ( grid_wl[:-1] + grid_wl[1:] ), grid_wl[-1] ]
    matrixW_dwl     = np.diff( matrixW_central )
    matrixW_sp_new    = np.dot( matrixW_convfunc, matrixW_sp * matrixW_dwl ) 

    return matrixW_sp_new


