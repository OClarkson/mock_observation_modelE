import numpy as np
import util_interp
import sys
import copy
from scipy import optimize

LOOPMAX = 100000


#=============================================================================
def RK4( list_dXds, array_X, ds ): 

    # X : l(s), z(h), beta, phi, xi
    n_var   = len( array_X )
    array_k = np.zeros([ 4, n_var ])

    # k0
    for ii in xrange( n_var ):
        array_k[0,ii] = list_dXds[ii]( array_X )
    # k1
    array_X_1 = array_X + 0.5*array_k[0]*ds
    for ii in xrange( n_var ):
        array_k[1,ii] = list_dXds[ii]( array_X_1 )
    # k2
    array_X_2 = array_X + 0.5*array_k[1]*ds
    for ii in xrange( n_var ):
        array_k[2,ii] = list_dXds[ii]( array_X_2 )
    # k3
    array_X_3 = array_X + array_k[2]*ds
    for ii in xrange( n_var ):
        array_k[3,ii] = list_dXds[ii]( array_X_3 )

    # step forward
    array_dX = ( array_k[0] + 2.*array_k[1] + 2.*array_k[2] + array_k[3] ) * ds / 6.

    return array_X + array_dX


#=============================================================================

#=============================================================================
def evolve_l_refraction( array_X, nn, dndr, r_planet, d_l ):

    # van der Werf (2008)
    # integration variable = s ( = l here )

    ll, zz, beta, phi = array_X
    inv_curvature = np.cos( beta ) / ( 1. + nn ) * dndr

    def ds_ds( array_X_tmp ):
        return 1.

    def dz_ds( array_X_tmp ):
        ll, zz, beta, phi = array_X_tmp
        return -1.*np.sin( beta )

    def dbeta_ds( array_X_tmp ):
        ll, zz, beta, phi = array_X_tmp
        return inv_curvature - 1. * np.cos( beta ) / ( r_planet + zz ) 

    def dphi_ds( array_X_tmp ):
        ll, zz, beta, phi = array_X_tmp
        return np.cos( beta ) / ( r_planet + zz )

    list_dXds = [ ds_ds, dz_ds, dbeta_ds, dphi_ds ]
    array_X_new   = RK4( list_dXds, array_X, d_l )

    return array_X_new


#=============================================================================
def refraction( layer_b, r_planet, dict_atmprof_funcZ, d_l ):

    b_num = len( layer_b )
    b_bottom_indx = -1

    z_top = layer_b[-1]

    matrixB_z_min   = np.zeros( b_num )
    matrixB_deflect = np.zeros( b_num )

    #------------------------------------------------
    # ray tracing
    #------------------------------------------------
    listB_func_lofz = []

    for ii in xrange( b_num - 1 ) :

        #------------------------------------------------
        # initialization
        # the planet is on the right-hand side
        # ray is traced from the planet
        #------------------------------------------------
        ll   = 0.

        # height
        zz   = z_top      

        # angle from the horizontal line to the ray
        z_i  = layer_b[ii]
        beta = np.arccos( ( r_planet + z_i ) / ( r_planet + z_top )  ) 

        # planet-centric angle
        phi  = -1. * beta 

        array_X = np.array([ ll, zz, beta, phi ])

        list_l = []
        list_z = []

        for loop in xrange( LOOPMAX ) :

            nn   = dict_atmprof_funcZ['refractivity']( array_X[1] )
            dndr = dict_atmprof_funcZ['dndr']( array_X[1] )

            array_X_new = evolve_l_refraction( array_X, nn, dndr, r_planet, d_l )
            zz_new      = array_X_new[1]
            beta_new    = array_X_new[2]

            ll, zz, beta, phi = array_X

            if ( beta_new < 0. ):

                # minimum altitude
                inv_curvature = np.cos( beta ) / ( 1. + nn ) * dndr
                d_l_adjust = beta / ( np.cos( beta ) / ( r_planet + zz ) + inv_curvature )
                ll, zz, beta, phi = evolve_l_refraction( array_X, nn, dndr, r_planet, d_l_adjust )
                list_l.append( ll )
                list_z.append( zz )
                l_max = ll

                matrixB_z_min[ii]   = zz
                matrixB_deflect[ii] = 2. * phi

                # l as a function of z
                array_l = l_max - np.array( list_l )
                array_z = np.array( list_z )

                func_lofz = util_interp.interp_1d_boundary( array_z[::-1], array_l[::-1], logx=False, logy=False, ext=3 )

                break

            if ( zz_new < 0. ):

                b_bottom_indx = ii
                func_lofz = 0.
                matrixB_deflect[ii] = -1

                break

            list_l.append( ll )
            list_z.append( zz )
            array_X = copy.deepcopy( array_X_new )

        if loop == LOOPMAX - 1 :
            util_errors.exit_msg('The ray is not fully traced. Increase LOOPMAX in ray_trace.py')

        listB_func_lofz.append( func_lofz )

    b_bottom_indx += 1

    return b_bottom_indx, matrixB_z_min, matrixB_deflect, listB_func_lofz

