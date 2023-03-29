#=============================================================================
# Module
#=============================================================================

import numpy as np
from scipy import integrate
from scipy  import constants
import util_interp
import util_errors
import cgs

#=============================================================================
# functions
#=============================================================================

#=============================================================================
def profile_nXSofZ_cld(layer_z, grid_wn, tuple_func_atmprof, Qext, alpha, fac):

    #------------------------------------------------
    # read water content profile
    #------------------------------------------------
    if (CLD_GCM):
        print('reading gcm output...')
        cld_profile = read_wcfile_gcm(FILE_HEIGHT, FILE_CLD, tuple_func_atmprof, ilat=23)
    else:
        cld_profile = read_wcfile(infile_cld)

    #------------------------------------------------
    # compute Qext, Qsca
    #------------------------------------------------
    cldlayer_nXSext = np.zeros(len(cld_profile['z']))
    cldlayer_nXSsca = np.zeros(len(cld_profile['z']))

    lattice_r = np.logspace(np.log10(D_MIN), np.log10(D_MAX), D_NUM)
    mesh_r,    mesh_zcld = np.meshgrid(lattice_r, cld_profile['z'])
    mesh_r,    mesh_lwc  = np.meshgrid(lattice_r, cld_profile['LWC'])
    mesh_dNdlogD = func_dNdlogD_gamma(mesh_r, D_EFF, mesh_lwc, alpha)
    mesh_nXSext  = np.pi*mesh_r**2*mesh_dNdlogD/mesh_r # integration with r
    cldlayer_nXSext = fac*Qext*integrate.trapz(mesh_nXSext, x=mesh_r, axis=1)

    func_kext = util_interp.interp_1d(cld_profile['z'], cldlayer_nXSext, logx=False, logy=False, bounds_error=False, fill_value=0.0)
#    func_ksca = util_interp.interp_1d(cld_profile['z'], nXSsca_cldlayers, order=1, logx=False, logy=False, fill_value=0.0)

    print("----------cloud data---------")
    for ii in range(len(cld_profile['z'])):
        print((cld_profile['z'][ii], cldlayer_nXSext[ii]))

    print("")
    layer_nXSext = func_kext(layer_z)

    print("")
    print("----------atmospheric data---------")
    for ii in range(len(layer_z)):
        print((layer_z[ii], layer_nXSext[ii]))


#    nXSsca_layers = func_ksca(grid_z)
    mesh_nXSext = np.tile(layer_nXSext, (len(grid_wn),1))
    return mesh_nXSext


#=============================================================================
#def read_wcfile(infile):
#    names = ['z', 'LWC', 'R_eff']
#    formats = ['f8', 'f8', 'f8']
#    data = np.loadtxt(infile, 
#                      dtype={'names':tuple(names),
#                             'formats':tuple(formats)}, 
#                      comments='#')
#    data['z']   = data['z']*1.e5 # km => cm
#    data['LWC']   = data['LWC']*1.e-6 # g/m3 => g/cm3
#    data['R_eff'] = data['R_eff']*1.e-4 # um => cm
#    return data

#=============================================================================
def read_wcfile_gcm(infile_height, infile_cld, tuple_func_atmprof, ilat=0):

    cld_profile = {}
    cld_profile['z']     = np.loadtxt(infile_height).T[ilat]*1.0e2
    cld_profile['LWC']   = np.loadtxt(infile_cld).T[ilat]

    print((cld_profile['z']))
    print((cld_profile['LWC']))

    #------------------------------------------------
    # check the order of altitude
    #------------------------------------------------
    reverse = 0
    z_old = cld_profile['z'][0]
    for zi in range(1,len(cld_profile['z'])):
        if (cld_profile['z'][zi] < z_old):
            reverse = reverse + 1
    if (reverse == len(cld_profile['z'])):
        cld_profile['z'] = cld_profile['z'][::-1]
        cld_profile['LWC'] = cld_profile['LWC'][::-1]
    elif (reverse != 0):
        errors.exit_msg("Check the order of FILE_CLD.")

    # actually, 'P', 'cldh2o', 'R_eff'
    func_TofZ, func_PofZ, func_MUofZ, dict_func_NofZ = tuple_func_atmprof
    cldlayer_P   = func_PofZ(cld_profile['z'])
    cldlayer_T   = func_TofZ(cld_profile['z'])
    cldlayer_MU  = func_MUofZ(cld_profile['z'])
    print(("cldlayer_P", cldlayer_P))
    print((cldlayer_P[0]/(cgs.RR*cldlayer_T[0])*constants.N_A))
    cldlayer_rho = cldlayer_MU*cldlayer_P/(cgs.RR*cldlayer_T)

    cld_profile['LWC'] = cld_profile['LWC']*cldlayer_rho

#    func_RHOofZ = util_interp.interp_1d_boundary(cld_profile['z'], cld_profile['LWC'], array_rho, logx=False, logy=False, order=1)
#

    return cld_profile



#=============================================================================
def func_NofD_exp(D, wc):

    # my thought
    # lambda0 = 
    # N0 = wc*lambda0**4.0/(np.pi*RHOL)
    # Kim & Del Genio (2013)
    N0      = 8.0 # 8.0 [1/cm4] = 8.e6 [1/m4]
    lambda0 = (np.pi*N0*RHOL/wc)**0.25
    NofD    = N0*np.exp(-lambda0*D)
    return NofD


#=============================================================================
def func_dNdlogD_gamma(D, r_eff, wc, alpha):

#    if (wc==0.0):
#        dNdlogD = 0.0*D
#    else:
    b = (alpha + 3.0)/r_eff
    factor = (4.0*np.pi*RHOL/3.0)
    def integrand(Rarg):
        return Rarg**alpha*np.exp(-b*Rarg)*Rarg**2
    a = wc/(factor*integrate.quad(integrand, D_MIN, D_MAX)[0])
    dNdlogD = a*D**alpha*np.exp(-b*D)
    return dNdlogD



