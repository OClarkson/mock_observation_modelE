import numpy as np
import sys
import cgs

from scipy import constants

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
params = {'mathtext.default': 'regular', 
          'font.size': 14 }
plt.rcParams.update(params)



R_PLANET = 6371. # km
R_STAR   = 0.3761 * 695700. # km
# R_STAR   = 0.2 * 695700. # km

H_MIN = 0.   # km
H_MAX = 63. # km

#---------------------------------
# plot
#---------------------------------

f, ax1 = plt.subplots(1, 1, figsize=(6,4))

ax1.set_ylim([H_MIN, H_MAX])
ax1.set_ylabel('effective altitude [km]')

ax1.set_xlabel(r'wavelength [$\mu $m]')
ax1.set_xlim( np.log10( [1, 5] ) )
plt.xticks( np.log10( np.array( [ 1, 2, 3, 4, 5 ] ) ), ['1', '2', '3', '4', '5' ] )
# ax1.set_xscale('log')

# list_color = [ '#000000', '#660000', '#CC0000', '#FF0000', '#FF9999', 'yellow' ]
list_color = [ '#660000', '#660000', '#FF0000', '#FF0000', '#FF9999', '#FF9999', 'yellow' ]
list_label = [ 'S$_X$=0.6', '', '  $_{ }$    1.0', '', ' $_{ }$     1.2', '' ]
list_linestyle = [ '-', '--', '-', '--', '-', '--' ]

#---------------------------------
# read .nc file
#---------------------------------
file_num = len( sys.argv ) - 1

for ii in xrange( file_num ):

    spfile_r = np.loadtxt( sys.argv[ii+1], unpack=True )
    ax1.plot( np.log10( spfile_r[0] ), spfile_r[1], linestyle=list_linestyle[ii], color=list_color[ii], label=list_label[ii] )

handles, labels = ax1.get_legend_handles_labels()
ax1.legend( handles, labels, fontsize=12, handletextpad=0.1, markerscale=0.5, loc=2 )


ax1.text( np.log10(3), 60, r'GJ 876 (M4V)', color='black', ha='left', va='top' )


ax2=ax1.twinx()

PPM_MIN = ( R_PLANET + H_MIN )**2 / ( R_STAR**2 ) * 1e6
PPM_MAX = ( R_PLANET + H_MAX )**2 / ( R_STAR**2 ) * 1e6
ax2.set_ylim([ PPM_MIN, PPM_MAX ])
ax2.set_ylabel('transit depth [ppm]')


plt.savefig( "temp.pdf", bbox_inches='tight')

