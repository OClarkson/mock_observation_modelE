#------------------------------------------------
# 1. input/output files
#------------------------------------------------

DATAFILE_DIR   = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/"
DATAFILE_TAG   = "1901.aijP1SoM40.nc"

RFILE          = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/P1SoM40.R"
SPECTRAL_FILE_SW   = 'spectral_files/sp_sw_ga7_dsa'
SPECTRAL_FILE_LW   = 'spectral_files/sp_lw_ga7_dsa'

OUTFILE_DIR    = 'out/'
OUTFILE_TAG    = "1901.aijP1SoM40_lw300_PHASEeq270"

#------------------------------------------------
# 2. geometry of observation
#------------------------------------------------

INC                =  90. # inclination of the planetary orbit with respect to line of sight [deg]
PHASE_EQ           = -90. # orbital phase of the spring equinox [deg]
SUBSTELLAR_LON_0   = -180. # longitude substellar point at the start of observation [deg]

FULL_PHASE         = True
if FULL_PHASE : 
    DIV_ORBIT      = 100.
else :
    PHASE0         = -90. # orbital phase of planet at the start of observation [deg]
    DT             =   1. # time interval of observation [hr]
    TIME_END       =  24. # duration of observation [hr]    



#------------------------------------------------
# 3. computation mode
#------------------------------------------------

MONTHLY = False

ShortWave_LightCurve = True # shortwave, lightcurve
ShortWave_Spectrum   = True # shortwave, spectrum
LongWave_LightCurve  = True # longwave,  lightcurve
LongWave_Spectrum    = True # longwave,  spectrum

#------------------------------------------------
# 4. Advanced
#------------------------------------------------

# Debugging mode ?
# If True, a new directory is not created
DEBUG_ON = False 





