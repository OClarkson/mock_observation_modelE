import netCDF4
import numpy as np
from copy import deepcopy

from molecules import *

import cgs
import constants
import util_errors
import sys

#=============================================================================
def extract_prof( infile, dict_NonCondensableGas, z_top, g_planet ) :

    print "Reading " + infile + " for atmospheric profile...   "

    # molecular weight of dry air
    mu_air_dry = 0.
    sum        = 0.
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
    for zi in xrange( 1, len( refrac_layers )-1 ):
        dndr_layers[zi] = ( refrac_layers[zi+1] - refrac_layers[zi-1] ) / ( dict_atmprof['z'][zi+1] - dict_atmprof['z'][zi-1] )
    dndr_layers[-1] = ( refrac_layers[-1] - refrac_layers[-2] ) / ( dict_atmprof['z'][-1] - dict_atmprof['z'][-2] )
    dict_atmprof['dndr'] = -1.*deepcopy( dndr_layers ) # flip the sign

    # polarizibility
    polar_layers = molecules['H2O']['polarizability'] * dict_atmprof['xH2O']
    sum_layers   = deepcopy( dict_atmprof['xH2O'] )
    for molename in dict_NonCondensableGas :
        if dict_NonCondensableGas[molename]=='otherwise' :
            polar_otherwise = molecules[molename]['polarizability']
        else :
            polar_layers += molecules[molename]['polarizability'] * ( 1. - dict_atmprof['xH2O'] ) * dict_NonCondensableGas[molename]
            sum_layers   += ( 1. - dict_atmprof['xH2O'] ) * dict_NonCondensableGas[molename]
    polar_layers += polar_otherwise * ( 1. - sum_layers )
    dict_atmprof['polarizability'] = polar_layers

    list_dict_atmprof = [ dict_atmprof ]
    list_theta        = [ 0. ]

            # print 'new'
            # for zi in xrange( len( dict_atmprof['z'] ) ):
            #     print theta/np.pi*180., dict_atmprof['z'][zi]*1e-5, dict_atmprof['temp'][zi]
            # print ''

    print "...Finish."


    return list_theta, list_dict_atmprof



