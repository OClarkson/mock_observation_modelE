import numpy as np
import util_interp

def read_O3file( filename, pres_levels ):

    p_layers, xO3_layers = np.loadtxt( filename, unpack=True )
    xO3_layers = xO3_layers * 1e-6 # ppm

    func_interp = util_interp.interp_1d( p_layers, xO3_layers )

    return func_interp( pres_levels )
