import numpy as np
import util_errors
import cgs
from molecules import *
import pandas as pd
import os
from scipy import interpolate


#=============================================================================
def nearestindex_WN(WN_lattice, wn) :
    """
    Returns the wavenumber index nearest to the given wavelength
    """
    idx = np.abs(WN_lattice-wn).argmin()
    return idx


#=============================================================================
def read_lookuptable( xsfile, grid_wn, grid_P, grid_T ):
    """
    Extract Look-up Table
    """

    print "Reading " + xsfile + "...   "


    if os.path.exists( xsfile ) : # check if the file exists

        data = np.load( xsfile )

        grid_P_org = data['P']*cgs.mbar_to_barye
        grid_T_org = data['T']
        grid_wn_org = data['WN']
        grid_XS_org = data['XS']

        if not ( np.allclose( grid_P_org, grid_P ) and np.allclose( grid_T_org, grid_T ) ) :
            util_errors.exit_msg("Grids on pressure or temperature are not consistent among cross section tables.")

        indx_min = nearestindex_WN( grid_wn_org, grid_wn[0]  )
        indx_max = nearestindex_WN( grid_wn_org, grid_wn[-1] )

        grid_XS = grid_XS_org[indx_min:indx_max+1]
        del grid_XS_org

        #### TEST (to be eventually removed)
        id1, id2, id3 = np.where(grid_XS <= 0)
        for i1, i2, i3 in zip(id1, id2, id3) :
            grid_XS[i1][i2][i3] = 1.e-48
        #### TEST

    else :

        grid_XS = np.zeros([ len( grid_wn ), len( grid_T ), len( grid_P ) ]) + 1.e-48

    return grid_XS


#=============================================================================
def griddata_line( list_mol, XSFILE_TAG, grid_wn, grid_T, grid_P, cnt_h2o_on=False ):
    
    # initialize the look-up table
    dict_griddata_logXSofWNTlogP = {}
    m_Tgrid, m_logPgrid = np.meshgrid( grid_T, np.log( grid_P ), indexing="ij" )
    dict_griddata_logXSofWNTlogP['coords'] = np.c_[ m_Tgrid.flatten(), m_logPgrid.flatten() ]

    # read lookup tables
    for molename in list_mol :

        if molecules[molename]['radiative'] :

            if molename == 'H2O' and cnt_h2o_on :
                xsfile = XSFILE_TAG + "H2O_c25.npz"
                XS_lookuptable     = read_lookuptable( xsfile, grid_wn, grid_P, grid_T )
                xsfile = XSFILE_TAG + "H2O_cnt.npz"
                XS_lookuptable_cnt = read_lookuptable( xsfile, grid_wn, grid_P, grid_T )
                dict_griddata_logXSofWNTlogP[molename] = np.log( XS_lookuptable + XS_lookuptable_cnt ).reshape([ len( grid_wn ), len( grid_T )*len( grid_P ) ])

                xsfile = XSFILE_TAG + "H2O-H2O_cnt.npz"
                XS_lookuptable_cnt2 = read_lookuptable( xsfile, grid_wn, grid_P, grid_T )
                dict_griddata_logXSofWNTlogP['H2O-H2O'] = np.log( XS_lookuptable_cnt2 ).reshape([ len( grid_wn ), len( grid_T )*len( grid_P ) ])

            else :
                xsfile = XSFILE_TAG + molename + ".npz"
                XS_lookuptable = read_lookuptable( xsfile, grid_wn, grid_P, grid_T )
                dict_griddata_logXSofWNTlogP[molename] = np.log( XS_lookuptable ).reshape([ len( grid_wn ), len( grid_T )*len( grid_P ) ])

    return dict_griddata_logXSofWNTlogP


#=============================================================================
def griddata_add_UV( molename, filename, grid_wn, grid_T, grid_P, dict_griddata_logXSofWNTP ):

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
    mesh_wn, mesh_t = np.meshgrid( grid_wn, grid_T, indexing='ij' )
    xi     = np.c_[ mesh_wn.flatten(), mesh_t.flatten() ]

    xi_logXS = interpolate.griddata( points, values, xi, method='nearest', fill_value=-48 )

    mesh_logXS = xi_logXS.reshape( [ len( grid_wn ), len( grid_T ), 1 ] )
    mesh_logXS = np.tile( mesh_logXS, len( grid_P ) ).reshape( [ len( grid_wn ), len( grid_T )*len( grid_P ) ] )
    tmp = np.exp( dict_griddata_logXSofWNTP[molename] ) + np.exp( mesh_logXS )
    dict_griddata_logXSofWNTP[molename] = np.log( tmp )

    return dict_griddata_logXSofWNTP





