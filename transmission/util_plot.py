import matplotlib.pyplot as plt



def myplot_func(filename, myfunc, x_min, x_max, xlog=False, ylog=False, str_title="", str_xlabel="x", str_ylabel="y", int_num=10) :
    """
    plot function
    """
    plt.title(str_title)
    plt.xlabel(str_xlabel, fontsize=20)
    plt.ylabel(str_ylabel, fontsize=20)
    plt.grid()
    x = np.linspace(x_min, x_max, num=int_num)
    if xlog : 
        x = np.logspace(np.log10(x_min), np.log10(x_max), num=int_num)
        plt.xscale('log')
    if ylog :
        plt.yscale('log')

    y = np.zeros(int_num)
    for i in range(int_num) :
        y[i] = myfunc(x[i])

    plt.tick_params(axis='both')
    plt.plot(x,y)
    plt.savefig(filename)


def myplot_tbl(filename, x_array, y_array, xlog=False, ylog=False, str_title="", str_xlabel="x", str_ylabel="y") :
    """
    plot tabulated data
    """

    plt.title(str_title)
    plt.xlabel(str_xlabel, fontsize=20)
    plt.ylabel(str_ylabel, fontsize=20)
    plt.grid()
    if xlog : 
        plt.xscale('log')
    if ylog :
        plt.yscale('log')

    plt.tick_params(axis='both')
    plt.plot(x_array,y_array)
    plt.savefig(filename)

