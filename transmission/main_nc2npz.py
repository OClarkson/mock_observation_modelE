import numpy as np
import netCDF4
import cgs
from scipy import interpolate

P_DIM = 54
T_DIM = 13

MU_H2O = 18. # g

XSFILE_CNTNM = "/Users/yuka/Dropbox/Abs_coeff/h2o_lbl_swf_pt702.nc"
XSFILE       = "xstbl/xstbl_HITRAN2012_00010-10000_m09991_H2O_c25.npz"
OUTFILE      = "xstbl/xstbl_HITRAN2012_00010-10000_m09991_H2O_cnt.npz"

#--------------------------------------------------------------------
def read_nc_cnt( infile ):

    ncfile_r = netCDF4.Dataset( infile,'r', format='NETCDF3_64BIT') 

    wn_grid = ncfile_r.variables['nu'][::100] 
    p_calc  = ncfile_r.variables['p_calc'][:]
    t_calc  = ncfile_r.variables['t_calc'][:]

    print 'shape of xs', ncfile_r.variables['kabs'][:,::100].T.shape

    xs_grid = ncfile_r.variables['kabs'][:,::100].T.reshape( [ len( wn_grid ), P_DIM, T_DIM ] )
    xs_grid = xs_grid.swapaxes(1,2)
    print 'xs_grid', xs_grid[0]

    ncfile_r.close()

    p_grid = p_calc[::T_DIM]
    t_grid = t_calc[:T_DIM]

    print 'p_grid', p_grid
    # unit conversion
    wn_grid = wn_grid * 1e-2  # m^-1 => cm^-1 
    p_calc  = p_calc  * 1e-2 # Pa => mbar
    xs_grid = xs_grid * 1e4 * 1e-3 # m2 kg-1 => cm2 / g^-1
    xs_grid = xs_grid * ( MU_H2O / cgs.NA )  # cm2 g-1 =>   cm2 / molecule
    print 'xs_grid2', xs_grid[0]

    return wn_grid, p_grid, t_grid, xs_grid


#--------------------------------------------------------------------
def interpolate_logXSofWNTP( infile, WN_lookuptable, TT_lookuptable, PP_lookuptable ) :

    wn_grid_cnt, p_grid_cnt, t_grid_cnt, xs_grid_cnt = read_nc_cnt( infile )

    print xs_grid_cnt.shape
    #### TEST (to be eventually removed)
    xs_grid_cnt[ np.where( xs_grid_cnt <= 1.e-50 ) ] = 1.e-50
    print 'xs_grid_cnt', xs_grid_cnt[0]
    #### TEST

    m_wn_cnt, m_T_cnt, m_logP_cnt = np.meshgrid( wn_grid_cnt, t_grid_cnt, np.log( p_grid_cnt ), indexing='ij' )
    coords = np.c_[ m_wn_cnt.flatten(), m_T_cnt.flatten(), m_logP_cnt.flatten() ]

    values = np.log( xs_grid_cnt.flatten() )

    m_wn, m_T, m_logP = np.meshgrid( WN_lookuptable, TT_lookuptable, np.log( PP_lookuptable ), indexing='ij' )
    flat_points = np.c_[ m_wn.flatten(), m_T.flatten(), m_logP.flatten() ]

#    flat_logXS_cnt  = interpolate.griddata( coords, values, flat_points, method='linear' ) 
    flat_logXS_cnt  = interpolate.griddata( coords, values, flat_points, method='nearest' ) 
    logXSofWNTlogP  = flat_logXS_cnt.reshape( [ len( WN_lookuptable ), len( TT_lookuptable ), len( PP_lookuptable ) ] )

    return logXSofWNTlogP



#=============================================================================
# main
#=============================================================================
if __name__ == "__main__":

    data = np.load( XSFILE )

    PP_grid     = data['P'] # mbar
    TT_grid     = data['T'] # K
    WN_grid     = data['WN'] # cm^-1
    XS_grid     = data['XS'] # cm^-2 

    logXSofWNTlogP = interpolate_logXSofWNTP( XSFILE_CNTNM, WN_grid, TT_grid, PP_grid )


    np.savez( OUTFILE, WN=WN_grid, T=TT_grid, P=PP_grid, XS=np.exp( logXSofWNTlogP ) )
