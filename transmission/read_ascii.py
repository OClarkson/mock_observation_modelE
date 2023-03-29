import netCDF4
import numpy as np
from copy import deepcopy

from molecules import *

import set_O3
from setup import l_O3, s_O3file

import cgs
import constants
import util_errors
import sys

#=============================================================================
def extract_prof( infile, dict_NonCondensableGas, p_min, p_max, g_planet ) :

    print(("Reading " + infile + " for atmospheric profile...   "))

    # molecular weight of dry air
    mu_air_dry = 0.
    sum        = 0.
    weight_otherwise = 0.
    for molename in dict_NonCondensableGas :
        if dict_NonCondensableGas[molename]=='otherwise' :
            weight_otherwise = molecules[molename]['weight']
        else :
            mu_air_dry += dict_NonCondensableGas[molename]*molecules[molename]['weight']
            sum        += dict_NonCondensableGas[molename]
    mu_air_dry += ( 1. - sum )*weight_otherwise

    # read ascii file
    atmprof = np.loadtxt( infile, unpack=True )

    dict_atmprof = {}
    dict_atmprof['z']    = atmprof[0][::-1]
    dict_atmprof['plm']  = atmprof[1][::-1]
    dict_atmprof['TempL'] = atmprof[2][::-1]
    dict_atmprof['q']    = atmprof[3][::-1]

    # O3
    if l_O3 :
        dict_atmprof['xO3']   = set_O3.read_O3file( s_O3file, dict_atmprof['plm']  )

    # unit conversion
    dict_atmprof['z']    = dict_atmprof['z']  *cgs.km_to_cm
    dict_atmprof['plm']  = dict_atmprof['plm']*cgs.mbar_to_barye

    # H2O vapor mixing ratio
    dict_atmprof['xH2O'] = dict_atmprof['q'] / constants.MU_H2O / ( dict_atmprof['q'] / constants.MU_H2O + ( 1 - dict_atmprof['q'] ) / mu_air_dry )

    # number density
    dict_atmprof['ndensity'] = dict_atmprof['plm'] / ( cgs.RR * dict_atmprof['TempL'] ) # mol

    # density
    layer_mu_atm = 1./ ( ( 1. - dict_atmprof['q'] ) / mu_air_dry + dict_atmprof['q'] / constants.MU_H2O )
    dict_atmprof['rho'] = dict_atmprof['ndensity'] * layer_mu_atm

    # refractivity
    refrac_layers  = molecules['H2O']['refractivity'] * ( dict_atmprof['ndensity'] * dict_atmprof['xH2O'] ) / cgs.amagat 

    sum_layers     = deepcopy( dict_atmprof['xH2O'] )
    refrac_otherwise = 0.
    for molename in dict_NonCondensableGas :
        if dict_NonCondensableGas[molename]=='otherwise' :
            refrac_otherwise = molecules[molename]['refractivity']
        else :
            refrac_layers += molecules[molename]['refractivity'] * ( dict_atmprof['ndensity'] * ( 1. - dict_atmprof['xH2O'] ) * dict_NonCondensableGas[molename] ) / cgs.amagat
            sum_layers    += ( 1. - dict_atmprof['xH2O'] ) * dict_NonCondensableGas[molename]
    refrac_layers += refrac_otherwise * ( dict_atmprof['ndensity'] * ( 1. - sum_layers ) ) / cgs.amagat
    dict_atmprof['refractivity'] = deepcopy( refrac_layers )

    # differential refractivity
    dndr_layers     = np.zeros_like( refrac_layers )
    dndr_layers[0]  = ( refrac_layers[1] - refrac_layers[0] ) / ( dict_atmprof['z'][1] - dict_atmprof['z'][0] )
    for zi in range( 1, len( refrac_layers )-1 ):
        dndr_layers[zi] = ( refrac_layers[zi+1] - refrac_layers[zi-1] ) / ( dict_atmprof['z'][zi+1] - dict_atmprof['z'][zi-1] )
    dndr_layers[-1] = ( refrac_layers[-1] - refrac_layers[-2] ) / ( dict_atmprof['z'][-1] - dict_atmprof['z'][-2] )
    dict_atmprof['dndr'] = -1.*deepcopy( dndr_layers ) # flip the sign

    list_dict_atmprof = [ dict_atmprof ]
    list_theta        = [ 0. ]

    return list_theta, list_dict_atmprof



