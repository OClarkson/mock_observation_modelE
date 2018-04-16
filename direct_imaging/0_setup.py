#------------------------------------------------
# 1. input/output files
#------------------------------------------------

s_RFile        = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/P1SoM40.R"
s_aijFile_Dir  = "/Users/yuka/Dropbox/GCM_RUNS/RUNS/P1SoM40/"

l_Monthly      = False
s_aijFile_Tag  = "1901.aijP1SoM40.nc" # if l_Monthly=True, specify the common suffix of the aij files, leaving out JAN/FEB/MAR etc. 
# s_aijFile_Tag  = "ANN1901.aijP1SoM40.nc" # if l_Monthly=False, specify the full file name

l_Socrates     = True
if l_Socrates :
    s_SpectralFile_SW = 'spectral_files/sp_sw_ga7_dsa'
    s_SpectralFile_LW = 'spectral_files/sp_lw_ga7_dsa'

s_outFile_Dir  = 'out/'
s_outFile_Tag  = "1901.aijP1SoM40_lw300_PHASEeq270"

#------------------------------------------------
# 2. geometry of observation
#------------------------------------------------

f_InclinationAngle_deg            =  90. # inclination of the planetary orbit with respect to line of sight [deg]
f_PhaseAngle_Equinox_deg          = -90. # orbital phase of the spring equinox [deg]
f_SubStellarLongitude_Initial_deg = -180. # longitude substellar point at the start of observation [deg]

#------------------------------------------------
# 3. computation mode
#------------------------------------------------

l_ShortWave_LightCurve = True # shortwave, lightcurve
l_ShortWave_Spectrum   = True # shortwave, spectrum
l_LongWave_LightCurve  = True # longwave,  lightcurve
l_LongWave_Spectrum    = True # longwave,  spectrum

l_FullOrbit = True
if l_FullOrbit : 
    i_DivideOrbit = 2000
else :
    f_PhaseAngle_Initial_deg = -90. # orbital phase of planet at the start of observation [deg]
    f_TimeInterval_hr        =   1. # time interval of observation [hr]
    f_TimeLimit_hr           =  24. # duration of observation [hr]    

#------------------------------------------------
# 4. Advanced
#------------------------------------------------

# Debugging mode ?
# If True, a new directory is not created
l_Debug = False 





