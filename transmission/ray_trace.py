import numpy as np
import util_interp
import sys

LOOPMAX = 100000


#=============================================================================
def evolve_l_refraction( ll, zz, beta, phi, nn, dndr, r_planet, d_l ):

    inv_curvature = np.cos( beta ) / ( 1. + nn ) * dndr

    d_z    = np.sin( beta ) * d_l
    d_phi  = np.cos( beta ) / ( r_planet + zz ) * d_l
    d_beta = inv_curvature * d_l - d_phi

    ll_new   = ll   + d_l
    zz_new   = zz   - d_z
    beta_new = beta + d_beta
    phi_new  = phi  + d_phi

    return ll_new, zz_new, beta_new, phi_new


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
        z_i  = layer_b[ii]

        # angle from the horizontal line to the ray
        beta = np.arccos( ( r_planet + z_i ) / ( r_planet + z_top )  ) 

        # planet-centric angle
        phi  = -1. * beta 

        # radius
        zz   = z_top      
        phi_ini = phi

        list_l = []
        list_z = []

        for loop in xrange( LOOPMAX ) :

            nn   = dict_atmprof_funcZ['refractivity']( zz )
            dndr = dict_atmprof_funcZ['dndr']( zz )
            ll_new, zz_new, beta_new, phi_new = evolve_l_refraction( ll, zz, beta, phi, nn, dndr, r_planet, d_l )

            if ( beta_new < 0. ):

                # minimum altitude
                inv_curvature = np.cos( beta ) / ( 1. + nn ) * dndr
                d_l_adjust = beta / ( np.cos( beta ) / ( r_planet + zz ) + inv_curvature )
                ll, zz, beta, phi = evolve_l_refraction( ll, zz, beta, phi, nn, dndr, r_planet, d_l_adjust )

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
            
            ll   = ll_new
            zz   = zz_new
            beta = beta_new
            phi  = phi_new

            list_l.append( ll )
            list_z.append( zz )

        if loop == LOOPMAX - 1 :
            util_errors.exit_msg('The ray is not fully traced. Increase LOOPMAX in ray_trace.py')

        listB_func_lofz.append( func_lofz )

    b_bottom_indx += 1

    return b_bottom_indx, matrixB_z_min, matrixB_deflect, listB_func_lofz

