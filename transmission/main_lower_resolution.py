#=============================================================================
#
# A code to lower the spectral resolution
# 
# Useage :
# python main_lower_resolution.py high-resolution-spectrum spectral-resolution
#
# Example :
# python main_lower_resolution.py out/test/mytest.txt 100.
#
#=============================================================================

from copy import deepcopy
import numpy as np
import sys

INFILE_TAG = sys.argv[1]
RESOL  = float(sys.argv[2])


#=============================================================================
def lower_resolution(wn_array, sp_array, resolution):

    count = 0

    wn_new_array = np.zeros(1)
    sp_new_array = np.zeros(1)

    wn_tmp = wn_array[0]
    jj = 0
    for ii in range(len(wn_array)):

        if wn_array[ii] > wn_tmp+wn_tmp/resolution :
            wn_new_array[jj] = wn_new_array[jj]/(count*1.0)
            sp_new_array[jj] = sp_new_array[jj]/(count*1.0)
            wn_new_array = np.r_[wn_new_array, np.zeros(1)]
            sp_new_array = np.r_[sp_new_array, np.zeros(1)]
            wn_tmp = wn_tmp + wn_tmp/resolution
            count = 0
            jj += 1

        wn_new_array[jj] += wn_array[ii]
        sp_new_array[jj] += sp_array[ii]
        count += 1

    wn_new_array[jj] = wn_new_array[jj]/(count*1.0)
    sp_new_array[jj] = sp_new_array[jj]/(count*1.0)

    return wn_new_array, sp_new_array



#=============================================================================
# main
#=============================================================================

if __name__ == "__main__":

    count = 0
    
    wn_new_array = np.zeros(1)
    sp_new_array = np.zeros(1)
        
    data = np.loadtxt( INFILE_TAG )
    wn_array = data.T[0]
    sp_array = data.T[1]
    if not ( 'wn' in sys.argv ) :
        wn_array = 1e4 / wn_array

    wn_new_array, sp_new_array = lower_resolution(wn_array, sp_array, RESOL)

    if not ( 'wn' in sys.argv ) :
        wn_new_array = 1e4 / wn_new_array
    
    data2 = np.dstack([wn_new_array, sp_new_array])
    np.savetxt(INFILE_TAG+"_R"+sys.argv[2], data2[0])

