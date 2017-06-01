import numpy as np
import errors
import netCDF4

deg2rad  = (np.pi/180.)
BAND_MAX = 100


#=======================================================================
def init_geometry( oblqty_deg, phase_eq_deg, phase0_deg, inc_deg ):

    oblqty   = deg2rad * oblqty_deg
    phase_eq = deg2rad * phase_eq_deg
    phase0   = deg2rad * phase0_deg
    inc      = deg2rad * inc_deg

    return oblqty, phase_eq, phase0, inc


#=======================================================================
def init_omega( p_spin, p_orbit ):

    omega_spin  = (2.0*np.pi/(p_spin))
    omega_orbit = (2.0*np.pi/(p_orbit))

    return omega_spin, omega_orbit



#=======================================================================
def init_area( nlat, nlon, lat_deg, lon_deg ):

    lat = lat_deg * deg2rad
    lon = lon_deg * deg2rad
    dlat = np.pi/nlat/2.0
    dlon = (2.*np.pi)/nlon
    area = ( np.sin( lat + dlat ) - np.sin( lat - dlat ) )*dlon
    return area


#=======================================================================
def get_weight( omega_spin, omega_orbit, oblqty, phase_eq, phase0, inc, 
                lat, lon, time ):
    """
    Calculate weight function
    """

    rot_x_oblqty   = np.array([[                     1.,                     0.,                 0. ],
                               [                     0.,       np.cos( oblqty ),   np.sin( oblqty ) ],
                               [                     0.,   -1.*np.sin( oblqty ),   np.cos( oblqty ) ]])

    rot_z_phase_eq = np.array([[     np.cos( phase_eq ),   -1.*np.sin( phase_eq ),                 0. ],
                               [     np.sin( phase_eq ),       np.cos( phase_eq ),                 0. ],
                               [                     0.,                       0.,                 1. ]])

    vecEO = np.array( [ np.sin( inc ), 0., np.cos( inc ) ] )
    lat    = deg2rad * lat
    lon    = deg2rad * lon + omega_spin * time

    vecERR = np.array([ np.cos( lat ) * np.cos( lon ),
                        np.cos( lat ) * np.sin( lon ),
                        np.sin( lat ) ])
    vecER  = np.dot( rot_z_phase_eq, np.dot( rot_x_oblqty, vecERR ) )

    phase  = phase0 + omega_orbit * time
    vecES  = np.array([ np.cos( phase ), np.sin( phase ), 0. ])
    cosTH0 = np.dot( vecES, vecER )
    cosTH1 = np.dot( vecEO, vecER )

    return cosTH0, cosTH1


