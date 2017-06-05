import cgs
import constants


#------------------------------------------------
# 1. input/output files
#------------------------------------------------

s_atmFile     = "data/us_standard.txt"

s_outFile_Dir = 'out/'
s_outFile_Tag = "us_standard"

#------------------------------------------------
# 2. look-up tables
#------------------------------------------------

s_xsFile_Tag  = "xstbl/xstbl_HITRAN2012_00010-10000_m09991_"

#------------------------------------------------
# 3. planetary system parameters (in cgs unit)
#------------------------------------------------

R_PLANET         = constants.R_earth
G_PLANET         = 9.8 * cgs.m_to_cm

R_STAR           = constants.R_sun
DISTANCE_TO_STAR = cgs.au_to_cm
IMPACT_PARAMETER = 0.0
ORBITAL_PHASE    = 0.0       # 0 at planetary transit.

DICT_NonCondensableGas = { 'O2' : 20.9476e-2, 
                           'CO2': 320.e-6, 
                           'CH4': 1.8e-6, 
                           'N2' : 'otherwise' }

#------------------------------------------------
# 4. computation mode
#------------------------------------------------

# Average atmosphere before calculation ?
l_atm_average    = False
if l_atm_average :
    i_atmave_num = 90

# Multicore ?
l_multicore      = False
if l_multicore :
    i_core_num   = 2

# Include Refraction ?
l_refraction     = True

# Include molecular absorption ?
l_molecular_absorption = True
if l_molecular_absorption :
    # H2O continuum ?
    l_H2O_continuum = True     

else :
    # Use user-defined wavenumber grids
    # (needed when gas absorption is not included
    f_Wavenumber_min = 1000
    f_Wavenumber_max = 10000
    i_Wavenumber_num = 100

# Include Rayleigh scattering ?
l_Rayleigh       = True

# Include cloud scattering ?
l_cloud          = False

if l_cloud :
    # Cloud particle grids
    f_cloud_Deff_liquid = 30.0e-4 # cm
    f_cloud_Deff_ice    = 30.0e-4 # cm
    f_factor            = 1.

# Make low-resolution spectra ?
l_lower_resolution  = True
if l_lower_resolution :
    f_resolution = 100

# Plot ?
l_Plot = True

#------------------------------------------------
# 5. computational parameters
# CHECK CONVERGENCE BY CHANGING THEM !!
#------------------------------------------------

# altitude of top of atmosphere
f_Z_top = 100.0e5 # cm

# grid point for integration
i_Z_num = 100

i_B_num = 100

# increment for ray tracing
f_dL    = 1.e5 # cm => 0.1 km

#------------------------------------------------
# 6. Advanced
#------------------------------------------------

# Debugging mode ?
# If True, a new directory is not created
l_Debug = False 



