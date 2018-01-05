import numpy as np
import util_errors
import cgs
from molecules import *
import pandas as pd
from scipy import interpolate


#=============================================================================
def nearestindex_WN(WN_lattice, wn) :
    """
    Returns the wavenumber index nearest to the given wavelength
    """
    idx = np.abs(WN_lattice-wn).argmin()
    return idx


#=============================================================================
def read_lookuptable( xsfile, wn_limit ):
    """
    Extract Look-up Table
    """

    print "Reading " + xsfile + "...   "

    data = np.load( xsfile )

    PP_grid     = data['P']*cgs.mbar_to_barye
    TT_grid     = data['T']
    WN_grid_org = data['WN']
    XS_grid_org = data['XS']

    #------------------------------------------------
    # check wavenumber range
    #------------------------------------------------
    if ( ( wn_limit[0] < min( WN_grid_org ) ) or ( wn_limit[1] > max( WN_grid_org ) ) ) : 
        util_errors.exit_msg("Wavenumbers set are out of range of look-up tables.")
    #------------------------------------------------

    indx_min = nearestindex_WN( WN_grid_org, wn_limit[0] )
    indx_max = nearestindex_WN( WN_grid_org, wn_limit[1] )

    WN_grid = WN_grid_org[indx_min:indx_max+1]
    del WN_grid_org
    XS_grid = XS_grid_org[indx_min:indx_max+1]
    del XS_grid_org

    #### TEST (to be eventually removed)
    id1, id2, id3 = np.where(XS_grid <= 0)
    for i1, i2, i3 in zip(id1, id2, id3) :
        XS_grid[i1][i2][i3] = 1.e-48
    #### TEST

    return WN_grid, TT_grid, PP_grid, XS_grid


#=============================================================================
def griddata_line( list_mol, XSFILE_TAG, wn_range, cnt_h2o_on=False ):

    wn_min, wn_max = wn_range

    # initialize cross section dictionary
    dict_griddata_logXSofWNTlogP = {}

    # read lookup tables
    for molename in list_mol :

        if molecules[molename]['radiative'] :

            if molename == 'H2O' and cnt_h2o_on :
                xsfile = XSFILE_TAG + "H2O_c25.npz"
                WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_lookuptable     = read_lookuptable(xsfile, (wn_min, wn_max))
                xsfile = XSFILE_TAG + "H2O_cnt.npz"
                WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_lookuptable_cnt = read_lookuptable(xsfile, (wn_min, wn_max))
                dict_griddata_logXSofWNTlogP[molename] = np.log( XS_lookuptable + XS_lookuptable_cnt ).reshape([ len( WN_lookuptable ), len( TT_lookuptable )*len( PP_lookuptable ) ])

                xsfile = XSFILE_TAG + "H2O-H2O_cnt.npz"
                WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_lookuptable_cnt2 = read_lookuptable(xsfile, (wn_min, wn_max))
                dict_griddata_logXSofWNTlogP['H2O-H2O'] = np.log( XS_lookuptable_cnt2 ).reshape([ len( WN_lookuptable ), len( TT_lookuptable )*len( PP_lookuptable ) ])

            else :
                xsfile = XSFILE_TAG + molename + ".npz"
                WN_lookuptable, TT_lookuptable, PP_lookuptable, XS_lookuptable = read_lookuptable(xsfile, (wn_min, wn_max))
                dict_griddata_logXSofWNTlogP[molename] = np.log( XS_lookuptable ).reshape([ len( WN_lookuptable ), len( TT_lookuptable )*len( PP_lookuptable ) ])


    m_Tgrid, m_logPgrid = np.meshgrid( TT_lookuptable, np.log(PP_lookuptable), indexing="ij" )
#    dict_griddata_logXSofWNTP['coords'] = np.c_[ m_WNgrid.flatten(), m_Tgrid.flatten(), m_logPgrid.flatten() ]
    dict_griddata_logXSofWNTlogP['coords'] = np.c_[ m_Tgrid.flatten(), m_logPgrid.flatten() ]

    return WN_lookuptable, dict_griddata_logXSofWNTlogP


#=============================================================================
def griddata_add_UV( molename, filename, grid_wn, dict_griddata_logXSofWNTP ):

    # read UV continuum absorption data
    print "Reading " + filename + " for " + molename + "...   "
    UVdata   = pd.read_table( filename, comment='#' )
    table_XS = UVdata.iloc[:,1:].as_matrix()
    table_T  = UVdata.columns[1:].astype( float ) 
    table_WL = UVdata.iloc[:,0].as_matrix()
    table_WN = 1e7 / table_WL

    # reverse 
    table_WN = table_WN[::-1]
    table_XS = table_XS[::-1,:]
    table_logXS = np.log( table_XS )

    mesh_WN, mesh_T = np.meshgrid( table_WN, table_T, indexing='ij' )
    points = np.c_[ mesh_WN.flatten(), mesh_T.flatten() ]
    values = table_logXS.flatten()

    # interpolate to match the grids of line absorption cross section    
    # -------------------------------------------
    # NOTE:
    # dimension of dict_interpolate.griddata_logXSofWNTP['XS']
    #  WN x T x P
    # -------------------------------------------
    grid_t = np.unique( dict_griddata_logXSofWNTP['coords'][:,0] )
    grid_p = np.unique( dict_griddata_logXSofWNTP['coords'][:,1] )
    mesh_wn, mesh_t = np.meshgrid( grid_wn, grid_t, indexing='ij' )
    xi     = np.c_[ mesh_wn.flatten(), mesh_t.flatten() ]

    xi_logXS = interpolate.griddata( points, values, xi, method='nearest', fill_value=-48 )

    mesh_logXS = xi_logXS.reshape( [ len( grid_wn ), len( grid_t ), 1 ] )
    mesh_logXS = np.tile( mesh_logXS, len( grid_p ) ).reshape( [ len( grid_wn ), len( grid_t )*len( grid_p ) ] )
    tmp = np.exp( dict_griddata_logXSofWNTP[molename] ) + np.exp( mesh_logXS )
    dict_griddata_logXSofWNTP[molename] = np.log( tmp )

    return dict_griddata_logXSofWNTP





