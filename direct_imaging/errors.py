import sys

def exit_msg(string):
    print("-----------------------------------------------------------------")
    print(" ERROR!!")
    print(" " + string)
    print("-----------------------------------------------------------------")
#    print "-------------------------- TERMINATION --------------------------"
    sys.exit()

def exit_longmsg(arraylike_string):
    print("-----------------------------------------------------------------")
    print(" ERROR!!")
    for ii in range(len(arraylike_string)):
        print(" " + arraylike_string[ii])
    print("-----------------------------------------------------------------")
#    print "-------------------------- TERMINATION --------------------------"
    sys.exit()

