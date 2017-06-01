#------------------------------------------------
# 1. input/output files
#------------------------------------------------

# DATAFILE_DIR   = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/"
# DATAFILE_TAG   = "1901.aijP1SoM40.nc"
DATAFILE_DIR       = "/Users/yuka/Dropbox/GCM_RUNS/ProxCenb_mway/"
DATAFILE_TAG       = "4000.aijProxCenb04b_TL.nc"

RFILE              = "/Users/yuka/Dropbox/GCM_RUNS/ProxCenb_mway/ProxCenb04b_TL.R"
# SPECTRAL_FILE  = '/Users/yuka/Dropbox/GCM_RUNS/RUNS/spectral_files/sp_lw_300_jm2'
SPECTRAL_FILE      = 'spectral_files/sp_sw_ga7_dsa'
SPECTRAL_FILE_LW   = 'spectral_files/sp_lw_ga7_dsa'

OUTFILE_DIR        = 'out/'
# OUTFILE_TAG    = "1901.aijP1SoM40_lw300_PHASEeq270_sw"
OUTFILE_TAG        = "4000.aijProxCenb04b_TL"

#------------------------------------------------
# 2. geometry of observation
#------------------------------------------------

INC                =   90. # inclination of the planetary orbit with respect to line of sight [deg]
PHASE_EQ           =  270. # orbital phase of the spring equinox [deg]
SUBSTELLAR_LON_0   = -180. # longitude substellar point at the start of observation [deg]

FULL_PHASE = True
if FULL_PHASE : 
    DIV_ORBIT      = 100.
else :
    DT             =   1. # time interval of observation [hr]
    PHASE0         = -90. # orbital phase of planet at the start of observation [deg]
    TIME_END       =  24. # duration of observation [hr]    

#------------------------------------------------
# 3. computation mode
#------------------------------------------------

# Shortwave (sw) or longwave (lw)
MODE = "sw"
# MODE = "lw"

MONTHLY = False

# Plot Spectrum and/or lightcurve
PLOT_LC = True
PLOT_SP = True

#------------------------------------------------
# 4. Advanced
#------------------------------------------------

# Debugging mode ?
# If True, a new directory is not created
DEBUG_ON = False 





