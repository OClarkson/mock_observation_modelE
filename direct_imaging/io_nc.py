import numpy as np
import errors
import netCDF4

deg2rad  = (np.pi/180.)
BAND_MAX = 300

#=======================================================================
def read_nc_socrates( infile, mode=False, param_in=False ):
    """
    To read netCDF file and extract albedo / outgoing flux
    """

    print "infile", infile

    # read file
    ncfile_r = netCDF4.Dataset( infile, 'r', format='NETCDF3_64BIT')        
    lat      = ncfile_r.variables['lat'][:]
    lon      = ncfile_r.variables['lon'][:]
    nlat     = len(lat)
    nlon     = len(lon)

    # mode
    if param_in :
        param = param_in
    elif ( mode == "SW" or mode == "sw" ):
        param = 'srup_toa_band_'
    elif ( mode == 'LW' or mode == "lw" ) :
        param = 'trup_toa_band_'
    else:
        errors.exit_msg( "Invalid mode ( SW of LW )" )

    # count number of bands
    for ii in xrange( 1, 1000 )  :
        if param+str(ii) not in ncfile_r.variables :
            band_num = ii
            break

    #--------------------------------
    lat2 = np.zeros_like(lat)
    for ilat in xrange( nlat ):
        lat2[ilat] = - 90.0 + ( 180.0 / nlat ) * ( ilat + 0.5 )
    lat = lat2
    #--------------------------------

    lon_mesh, lat_mesh = np.meshgrid( lon, lat )
    lat_flatten = lat_mesh.flatten()
    lon_flatten = lon_mesh.flatten()

    array_data   = np.zeros( [ len( lat_flatten ), band_num - 1 ] )

    for ii in xrange( 1, band_num ):

        array_data[:,ii-1] = ncfile_r.variables[param+str(ii)][:].flatten()

        if ( mode == "SW" or mode == "sw" ):

            array_denomi = ncfile_r.variables['srdn_toa_band_'+str(ii)][:].flatten()

            # when incident radiation is zero, it is set to 1 and albedo is zero. 
            if len( np.where( array_denomi <= 0.)[0] ) > 0 :
                array_denomi[np.where( array_denomi <= 0. )] = 1.
                array_data[np.where( array_denomi <= 0. ),ii-1] = 0.

            array_data[:,ii-1] = array_data[:,ii-1] / array_denomi

    return nlat, nlon, lat_flatten, lon_flatten, array_data


#=======================================================================
def read_nc_giss( infile, mode=False, param_in=False ):
    """
    To read netCDF file and extract albedo / outgoing flux
    """

    print "infile", infile

    # read file
    ncfile_r = netCDF4.Dataset( infile, 'r', format='NETCDF3_64BIT')
    lat      = ncfile_r.variables['lat'][:]
    lon      = ncfile_r.variables['lon'][:]
    nlat     = len(lat)
    nlon     = len(lon)

    # mode
    if param_in :
        param = param_in
    elif ( mode == "SW" or mode == "sw" ):
        param = 'plan_alb'
    elif ( mode == 'LW' or mode == "lw" ) :
        param = 'trnf_toa'
    else:
        errors.exit_msg( "Invalid mode ( SW of LW )" )

    #--------------------------------
    lat2 = np.zeros_like(lat)
    for ilat in xrange( nlat ):
        lat2[ilat] = - 90.0 + ( 180.0 / nlat ) * ( ilat + 0.5 )
    lat = lat2
    #--------------------------------

    lon_mesh, lat_mesh = np.meshgrid( lon, lat )
    lat_flatten = lat_mesh.flatten()
    lon_flatten = lon_mesh.flatten()


    if ( mode == "SW" or mode == "sw" ):
        array_data_tmp = ncfile_r.variables[param][:]
        array_data_tmp = np.maximum( array_data_tmp, 0. ) * 0.01

    elif ( mode == 'LW' or mode == "lw" ) :
        array_data_tmp = ncfile_r.variables[param][:] * -1.

    array_data      = np.zeros( [ len( lat_flatten ), 1 ] )
    array_data[:,0] = array_data_tmp.flatten()

    return nlat, nlon, lat_flatten, lon_flatten, array_data



