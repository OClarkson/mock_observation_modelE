import cgs
import constants


#------------------------------------------------
# 1. input/output files
#------------------------------------------------

DATAFILE         = "data/us_standard.txt"
OUTFILE_DIR      = "out/"
OUTFILE_TAG      = "us_standard"

#------------------------------------------------
# 2. look-up tables
#------------------------------------------------

XSFILE_TAG       = "xstbl/xstbl_HITRAN2012_00010-10000_m09991_"

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
ATMAVE_ON        = False
if ATMAVE_ON :
    ATMAVE_NUM   = 90

# Multicore ?
MULTICORE_ON     = False
if MULTICORE_ON :
    CORE_NUM     = 2

# Include Refraction ?
REFRACTION_ON    = True

# Include molecular absorption ?
MOLABS_ON        = True
if MOLABS_ON :
    # H2O continuum ?
    CNTNM_H2O_ON = True     

else :
    # Use user-defined wavenumber grids
    # (needed when gas absorption is not included
    WN_MIN       = 1000
    WN_MAX       = 10000
    WN_NUM       = 100

# Include Rayleigh scattering ?
RAYLEIGH_ON      = True

# Include cloud scattering ?
CLD_ON           = False

# Cloud particle grids
CLD_D_EFF_LIQUID = 30.0e-4 # cm
CLD_D_EFF_ICE    = 30.0e-4 # cm
FACTOR           = 1.

# Make low-resolution spectra ?
LOW_RES_ON       = True
if LOW_RES_ON :
    RESOLUTION   = 100


#------------------------------------------------
# 5. computational parameters
# CHECK CONVERGENCE BY CHANGING THEM !!
#------------------------------------------------

# altitude of top of atmosphere
Z_TOP            = 100.0e5 # cm

# grid point for integration
Z_NUM            = 100

B_NUM            = 100

# increment for ray tracing
D_L              = 1.e5 # cm => 0.1 km


#------------------------------------------------
# 6. Advanced
#------------------------------------------------

# Debugging mode ?
# If True, a new directory is not created
DEBUG_ON = False 



