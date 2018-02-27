import numpy  as np
from scipy import interpolate
import cgs
import sys
from copy import deepcopy

from constants import RHO_H2O

from molecules import molecules

from setup import l_cloud, l_O3
if l_cloud :
    from setup import f_cloud_Deff_liquid, f_cloud_Deff_ice, f_factor

#=============================================================================
def get_nXS_molabs( layer_z, grid_wn, dict_griddata_logXSofWNTP, dict_NonCondensableGas, dict_atmprof_funcZ ) :

    matrixZ_pres = dict_atmprof_funcZ['plm'  ]( layer_z )
    matrixZ_temp = dict_atmprof_funcZ['TempL']( layer_z )
    matrixZ_xH2O = dict_atmprof_funcZ['xH2O' ]( layer_z )
    if l_O3 :
        matrixZ_xO3  = dict_atmprof_funcZ['xO3'  ]( layer_z )
    else :
        matrixZ_xO3 = np.zeros_like( matrixZ_xH2O )
    flat_points = np.c_[ matrixZ_temp, np.log( matrixZ_pres ) ]

    matrixWZ_nXS = np.zeros( [ len( grid_wn ), len( layer_z ) ] )

    list_mol = dict_griddata_logXSofWNTP.keys()
    list_mol.remove('coords')

    #---------- map function -----------
    def expand_wn( ii ):

        flat_logXS  = interpolate.griddata( dict_griddata_logXSofWNTP['coords'], dict_griddata_logXSofWNTP[molename][ii], flat_points, method='linear' ) 
        flat_XS = np.exp( flat_logXS )
        flat_XS[np.isnan( flat_XS )] = 0.
        return flat_XS

    #-----------------------------------

    for molename in list_mol :

        matrixWZ_XS = map( expand_wn, np.arange( len( grid_wn ) ) )
        matrixWZ_XS = np.array( matrixWZ_XS )

        if molename == 'H2O-H2O' :
            matrixZ_rho         = dict_atmprof_funcZ['ndensity']( layer_z ) * matrixZ_xH2O * molecules['H2O']['weight']
            matrixWZ_nXS_each = matrixWZ_XS * matrixZ_rho**2
            matrixWZ_nXS     += matrixWZ_nXS_each

        else :
            if molename == 'H2O' :
                matrixZ_n = dict_atmprof_funcZ['ndensity']( layer_z ) * cgs.NA * matrixZ_xH2O
            elif molename == 'O3' :
                matrixZ_n = dict_atmprof_funcZ['ndensity']( layer_z ) * cgs.NA * matrixZ_xO3
            else :
                matrixZ_n = dict_atmprof_funcZ['ndensity']( layer_z ) * cgs.NA * ( 1. - matrixZ_xH2O - matrixZ_xO3 ) * dict_NonCondensableGas[molename]

            matrixWZ_nXS_each = matrixWZ_XS * matrixZ_n 
            matrixWZ_nXS     += matrixWZ_nXS_each

    return matrixWZ_nXS


#=============================================================================
def get_nXS_Rayleigh( layer_z, grid_wn, dict_NonCondensableGas, dict_atmprof_funcZ ) :
    """
    Rayleigh scattering
    """

    matrixWZ_wn, matrixWZ_z = np.meshgrid( grid_wn, layer_z, indexing='ij' )
    matrixWZ_wl   = 1.0 / matrixWZ_wn
    matrixWZ_nXS  = np.zeros_like( matrixWZ_wn )

    matrixZ_pres  = dict_atmprof_funcZ['plm'     ]( layer_z )
    matrixZ_temp  = dict_atmprof_funcZ['TempL'   ]( layer_z )
    matrixZ_xH2O  = dict_atmprof_funcZ['xH2O'    ]( layer_z )

    if l_O3 :
        matrixZ_xO3   = dict_atmprof_funcZ['xO3'     ]( layer_z )
    else :
        matrixZ_xO3 = np.zeros_like( matrixZ_xH2O )

    matrixZ_n0    = dict_atmprof_funcZ['ndensity']( layer_z ) * cgs.NA
    matrixZ_alpha = np.zeros_like( layer_z )

    if l_O3 :
        list_mol = dict_NonCondensableGas.keys() + ['H2O', 'O3']
    else :
        list_mol = dict_NonCondensableGas.keys() + ['H2O']

    for molename in list_mol :

        if molename == 'H2O' :
            matrixZ_n = matrixZ_n0 * matrixZ_xH2O
        elif molename == 'O3' :
            matrixZ_n = matrixZ_n0 * matrixZ_xO3
        elif dict_NonCondensableGas[molename] == 'otherwise' :
            tmp_NonCondensableGas = deepcopy( dict_NonCondensableGas )
            tmp_NonCondensableGas[molename] = 0.
            sum_mixingratio = sum( tmp_NonCondensableGas.values() )
            matrixZ_n = matrixZ_n0 * ( 1. - matrixZ_xH2O - matrixZ_xO3 ) * ( 1. - sum_mixingratio )
        else :
            matrixZ_n = matrixZ_n0 * ( 1. - matrixZ_xH2O - matrixZ_xO3 ) * dict_NonCondensableGas[molename]
                    
        matrixZ_refract = molecules[molename]['refractivity'] * ( matrixZ_n / cgs.NA ) / cgs.amagat

        nonzero_indx = np.where( matrixZ_n > 0. )[0]
        matrixZ_alpha = ( 3. / ( 4. * np.pi * matrixZ_n[nonzero_indx] )) * ( ( matrixZ_refract[nonzero_indx] + 1. )**2 -1 ) / ( ( matrixZ_refract[nonzero_indx] + 1. )**2 + 2 )
        matrixWZ_nXS[:,nonzero_indx]  += 128.0 * np.pi**5 / ( 3.0 * matrixWZ_wl[:,nonzero_indx]**4 ) * matrixZ_n[nonzero_indx] * matrixZ_alpha**2

    return matrixWZ_nXS


#=============================================================================
def get_nXS_cld( layer_z, grid_wn, dict_atmprof_funcZ ) :
    """
    Opacity due to cloud particles
    """


    CCC = 2.* f_factor

    matrixWZ_wn, matrixWZ_z = np.meshgrid( grid_wn, layer_z, indexing='ij' )

    matrixWZ_rho_atm = dict_atmprof_funcZ['rho']( matrixWZ_z ) 

    # ice
    matrixWZ_icecld      = dict_atmprof_funcZ['icecld']( matrixWZ_z )
    matrixWZ_rho_cld_ice = matrixWZ_rho_atm * matrixWZ_icecld
    matrixWZ_nXS         = 3. * CCC * matrixWZ_rho_cld_ice / ( 4. * f_cloud_Deff_ice ) / RHO_H2O

    # liquid
    matrixWZ_wtrcld         = dict_atmprof_funcZ['wtrcld']( matrixWZ_z )
    matrixWZ_rho_cld_liquid = matrixWZ_rho_atm * matrixWZ_wtrcld
    matrixWZ_nXS           += 3. * CCC * matrixWZ_rho_cld_liquid / ( 4. * f_cloud_Deff_liquid ) / RHO_H2O

    return matrixWZ_nXS


