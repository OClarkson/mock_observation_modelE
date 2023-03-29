import numpy as np

#------------------------------------
# dynamical parameters
#------------------------------------
param_default = { 'obliquity': 23.44,                      # planetary obliquity [deg]
                  'siderealrotationperiod': 86636.7123288, # spin period [sec] 
                  'siderealorbitalperiod' : 31536000. }    # orbital period [sec]


#=======================================================================
def extract_param( file, list_param, type='txt' ):
    
    # R file
    ld = open( file )
    lines = ld.readlines()
    ld.close()

    list_value = []
    for param in list_param :

        # default
        if param in list(param_default.keys()):
            value = param_default[param]
        else:
            value = 0.

        # search in R
        for line0 in lines :
            line = line0.lower()
            if line.find( param.lower() ) >= 0 :
                if line[0] !='!' :
                    start = line.find( param.lower() ) + len( param ) + 1
                    end   = min( line.find(' '), line.find('\n'), line.find('!'))
                    if type=='float' :
                        print('line[start:end]', line[start:end])
                        str = line[start:end].replace( 'd', 'e' )
                        value = float( str )
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
    for ll in range( len( lines ) ):
        if lines[ll].find( start ) >= 0 :
            l_start = ll
            while lines[ll].find( end ) < 0 :
                ll = ll + 1
            l_end = ll

    data = np.genfromtxt( file, skip_header=l_start+4, skip_footer=len(lines)-l_end )

    return data

