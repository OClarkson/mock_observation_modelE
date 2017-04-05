#------------------------------------------------
# 1. input/output files
#------------------------------------------------

DATAFILE_DIR   = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/"
DATAFILE_TAG   = "1901.aijP1SoM40.nc"
RFILE          = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/P1SoM40.R"
SPECTRAL_FILE  = 'spectral_file/sp_lw_300_jm2'

OUTFILE_DIR    = 'out/'
OUTFILE_TAG    = "1901.aijP1SoM40_lw300_PHASEeq270"


#------------------------------------------------
# 2. geometry of observation
#------------------------------------------------

INC            =  90. # inclination of the planetary orbit with respect to line of sight [deg]
PHASE_EQ       = -90. # orbital phase of the spring equinox [deg]
PHASE0         = -90. # orbital phase of planet at the start of observation [deg]
TIME_END       =  24. # duration of observation [hr]
DT             =   1. # time interval of observation [hr]


#------------------------------------------------
# 3. computation mode
#------------------------------------------------

# Shortwave (sw) or longwave (lw)
# MODE = "sw"
MODE = "lw"

# Plot Spectrum and/or lightcurve
PLOT_LC = True
PLOT_SP = True


#------------------------------------------------
# 4. Advanced
#------------------------------------------------

# Debugging mode ?
# If True, a new directory is not created
DEBUG_ON = False 





