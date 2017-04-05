# wavelength :  589.3 nm
# 

# refractive index 
# http://www.kayelaby.npl.co.uk/general_physics/2_5/2_5_7.html

# polarizability
# http://cccbdb.nist.gov/polcalccomp2.asp?method=3&basis=4

molecules = {
    'N2'  : { 'weight'         : 28.014, 
              'refractivity'   : 0.000298, 
              'radiative'      : False,
              'polarizability' : 1.710e-24 }, # for now, N2-N2 molecular absorption is not included
    'O2'  : { 'weight'         : 32.00,
              'refractivity'   : 0.000271, 
              'radiative'      : True,
              'polarizability' : 1.562e-24 }, 
    'O3'  : { 'weight'         : 47.998,
              'refractivity'   : 0.,  # ???
              'radiative'      : True,
              'polarizability' : 3.079e-24 }, 
    'N2O' : { 'weight'         : 44.0128, 
              'refractivity'   : 0.000516, 
              'radiative'      : True,
              'polarizability' : 2.998e-24 }, 
    'CO2' : { 'weight'         : 44.01, 
              'refractivity'   : 0.000449,  
              'radiative'      : True,
              'polarizability' : 2.507e-24 }, 
    'H2O' : { 'weight'         : 18.01, 
              'refractivity'   : 0.000256, 
              'radiative'      : True,
              'polarizability' : 1.501e-24 },
    'CH4' : { 'weight'         : 16.042, 
              'refractivity'   : 0.000444, 
              'radiative'      : True,
              'polarizability' : 2.448e-24 },
    'CO'  : { 'weight'         : 28.01, 
              'refractivity'   : 0.000338, 
              'radiative'      : True,
              'polarizability' : 1.953e-24 }
    }
