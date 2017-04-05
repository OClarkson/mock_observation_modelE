import numpy as np


#=============================================================================
def make_wngrid( wn_min, wn_max, resolution ):

    #    R = ( wn + dwn / 2 ) /  dwn
    # => dwn = wn / ( R - 1./2. )
    list_wn_edge   = [ wn_min ]
    list_wn_center = [  ]
    wn = wn_min
    while wn < wn_max :
        dwn = wn / ( resolution - 0.5 )
        list_wn_center.append( wn + 0.5*dwn  )
        wn = wn + dwn
        list_wn_edge.append( wn )
    grid_wn_edge   = np.array( list_wn_edge   )
    grid_wn_center = np.array( list_wn_center )
    return grid_wn_edge, grid_wn_center


#=============================================================================
def lower_resolution( grid_wn, grid_sp, resolution ):

    grid2_wn_edge, grid2_wn_center = make_wngrid( grid_wn[0], grid_wn[-1], resolution )

    jj = 0
    sp_new_array = np.zeros( len( grid2_wn_edge ) - 1 )
    count = 0
    for ii in range(len(grid_wn)):

        if grid_wn[ii] > grid2_wn_edge[jj+1] :
            sp_new_array[jj] = sp_new_array[jj]/(count*1.0)
            jj    = np.where( grid_wn[ii] > grid2_wn_edge )[0][-1]
            count = 0

        sp_new_array[jj] += grid_sp[ii]
        count += 1

    sp_new_array[jj] = sp_new_array[jj]/(count*1.0)

    return grid2_wn_center, sp_new_array



