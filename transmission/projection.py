import numpy as np

def projected_on_stellardisk( theta, z_i, deflection_angle, dict_geom ) :

    X0 = dict_geom['distance_to_star'] * np.sin( dict_geom['orbital_phase'] )
    Y0 = dict_geom['r_star'] * dict_geom['impact_parameter']

    alpha = dict_geom['distance_to_star'] * np.tan( deflection_angle )
    dx    = -1. * ( dict_geom['r_planet'] + z_i - alpha ) * np.sin( theta )
    dy    =       ( dict_geom['r_planet'] + z_i - alpha ) * np.cos( theta )

    x_normalized = ( X0 + dx ) / dict_geom['r_star']
    y_normalized = ( Y0 + dy ) / dict_geom['r_star']

    return np.array( [ x_normalized, y_normalized ] )


def projected_on_planetdisk( theta, z_i, dict_geom ) :

    X0 = dict_geom['distance_to_star'] * np.sin( dict_geom['orbital_phase'] )
    Y0 = dict_geom['r_star'] * dict_geom['impact_parameter']

    dx    = -1. * ( dict_geom['r_planet'] + z_i ) * np.sin( theta )
    dy    =       ( dict_geom['r_planet'] + z_i ) * np.cos( theta )

    x_normalized = ( X0 + dx ) / dict_geom['r_star']
    y_normalized = ( Y0 + dy ) / dict_geom['r_star']

    return np.array( [ x_normalized, y_normalized ] )



def area_on_stellardisk( theta, d_theta, layer_z, layer_deflect, dict_geom ) :

    X0 = dict_geom['distance_to_star'] * np.sin( dict_geom['orbital_phase'] ) / dict_geom['r_star']
    Y0 = dict_geom['r_star'] * dict_geom['impact_parameter'] / dict_geom['r_star']

    layer_alpha = dict_geom['distance_to_star'] * np.tan( layer_deflect )
    layer_gamma = ( dict_geom['r_planet'] + layer_z - layer_alpha ) / dict_geom['r_star']
    layer_dx    = -1. * layer_gamma * np.sin( theta )
    layer_dy    =       layer_gamma * np.cos( theta )

    layer_gamma_edge = ( layer_gamma[1:] + layer_gamma[:-1] ) / 2.
    layer_gamma_edge = np.r_[ layer_gamma[0], layer_gamma_edge, layer_gamma[-1] ]
    layer_dgamma     = np.fabs( np.diff( layer_gamma_edge ) )

    layer_x = ( X0 + layer_dx ) 
    layer_y = ( Y0 + layer_dy )
    layer_d = layer_x**2 + layer_y**2

    layer_area = np.fabs( layer_gamma ) * layer_dgamma * d_theta
    layer_area[ np.where( layer_d > 1. ) ] = 0.

#    for zi in xrange( len( layer_gamma ) ):
#        print layer_z[zi], layer_dx[zi], layer_dy[zi], layer_deflect[zi], layer_gamma[zi], layer_area[zi]
#    print ''

    layer_area[ np.where( layer_deflect < 0. ) ] = 0.

#    for zi in xrange( len( layer_gamma ) ):
#        print layer_z[zi], layer_dx[zi], layer_dy[zi], layer_deflect[zi], layer_gamma[zi], layer_area[zi]
#    print ''


    return layer_area



def area_on_stellardisk_2( theta, d_theta, layer_z, layer_deflect, dict_geom ) :

    X0 = dict_geom['distance_to_star'] * np.sin( dict_geom['orbital_phase'] ) / dict_geom['r_star']
    Y0 = dict_geom['r_star'] * dict_geom['impact_parameter'] / dict_geom['r_star']

#    layer_alpha = dict_geom['distance_to_star'] * np.tan( layer_deflect )
#    layer_gamma = ( dict_geom['r_planet'] + layer_z - layer_alpha ) / dict_geom['r_star']
    layer_gamma = ( dict_geom['r_planet'] + layer_z ) / dict_geom['r_star']
    layer_dx    = -1. * layer_gamma * np.sin( theta )
    layer_dy    =       layer_gamma * np.cos( theta )

    layer_gamma_edge = ( layer_gamma[1:] + layer_gamma[:-1] ) / 2.
    layer_gamma_edge = np.r_[ layer_gamma[0], layer_gamma_edge, layer_gamma[-1] ]
    layer_dgamma     = np.fabs( np.diff( layer_gamma_edge ) )

    layer_x = ( X0 + layer_dx ) 
    layer_y = ( Y0 + layer_dy )
    layer_d = layer_x**2 + layer_y**2

    layer_area = np.fabs( layer_gamma ) * layer_dgamma * d_theta
    layer_area[ np.where( layer_d > 1. ) ] = 0.

#    for zi in xrange( len( layer_gamma ) ):
#        print layer_z[zi], layer_dx[zi], layer_dy[zi], layer_deflect[zi], layer_gamma[zi], layer_area[zi]
#    print ''

    layer_area[ np.where( layer_deflect < 0. ) ] = 0.

#    for zi in xrange( len( layer_gamma ) ):
#        print layer_z[zi], layer_dx[zi], layer_dy[zi], layer_deflect[zi], layer_gamma[zi], layer_area[zi]
#    print ''

    return layer_area


def area_on_stellardisk( theta, d_theta, layer_z, layer_deflect, dict_geom ) :

    X0 = dict_geom['distance_to_star'] * np.sin( dict_geom['orbital_phase'] ) / dict_geom['r_star']
    Y0 = dict_geom['impact_parameter']

    layer_alpha = dict_geom['distance_to_star'] * np.tan( layer_deflect )
    layer_gamma = ( dict_geom['r_planet'] + layer_z - layer_alpha ) / dict_geom['r_star']
    layer_dx    = -1. * layer_gamma * np.sin( theta )
    layer_dy    =       layer_gamma * np.cos( theta )

    layer_gamma2 = ( dict_geom['r_planet'] + layer_z ) / dict_geom['r_star']

    layer_gamma2_edge = ( layer_gamma2[1:] + layer_gamma2[:-1] ) / 2.
    layer_gamma2_edge = np.r_[ layer_gamma2[0], layer_gamma2_edge, layer_gamma2[-1] ]
    layer_dgamma2     = np.fabs( np.diff( layer_gamma2_edge ) )

    layer_x = ( X0 + layer_dx ) 
    layer_y = ( Y0 + layer_dy )
    layer_d = layer_x**2 + layer_y**2

#    for zi in xrange( len( layer_z ) ):
#        print layer_x[zi], layer_y[zi], layer_d[zi], layer_alpha[zi]/dict_geom['r_star']

    layer_area = np.fabs( layer_gamma2 ) * layer_dgamma2 * d_theta
    layer_area[ np.where( layer_d > 1. ) ] = 0.

#    for zi in xrange( len( layer_gamma ) ):
#        print layer_z[zi], layer_dx[zi], layer_dy[zi], layer_deflect[zi], layer_gamma[zi], layer_area[zi]
#    print ''

#    layer_area[ np.where( layer_deflect < 0. ) ] = 0.

#    for zi in xrange( len( layer_gamma ) ):
#        print layer_z[zi], layer_dx[zi], layer_dy[zi], layer_deflect[zi], layer_gamma[zi], layer_area[zi]
#    print ''

    return layer_area
