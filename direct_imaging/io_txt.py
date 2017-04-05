import numpy as np

#------------------------------------
# dynamical parameters
#------------------------------------
param_default = { 'OBLIQUITY': 23.44,                      # planetary obliquity [deg]
                  'siderealRotationPeriod': 86636.7123288, # spin period [sec] 
                  'siderealOrbitalPeriod' : 31536000. }    # orbital period [sec]


#=======================================================================
def extract_param( file, list_param, type='txt' ):
    
    # R file
    ld = open( file )
    lines = ld.readlines()
    ld.close()

    list_value = []
    for param in list_param :

        # default
        if param in param_default.keys():
            value = param_default[param]
        else:
            value = 0.

        # search in R
        for line in lines:
            if line.find( param ) >= 0 :
                if line[0] !='!' :
                    start = line.find( param ) + len( param ) + 1
                    end   = line.find(' ')
                    if type=='float' :
                        value = float( line[start:end] )
                    else :
                        value = line[start:end]

        list_value.append( value )

    return list_value


#=======================================================================
def extract_block( file, start, end ):

    # R file
    ld = open( file )
    lines = ld.readlines()
    ld.close()

    # search
    for ll in xrange( len( lines ) ):
        if lines[ll].find( start ) >= 0 :
            l_start = ll
            while lines[ll].find( end ) < 0 :
                ll = ll + 1
            l_end = ll

    data = np.genfromtxt( file, skip_header=l_start+1, skip_footer=len(lines)-l_end )

    return data

