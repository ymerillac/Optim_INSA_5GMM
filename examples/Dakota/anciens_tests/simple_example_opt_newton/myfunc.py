import sys
import numpy
import dakota_utils #--> le module a creer

def my_func(x):
    return x[0]**2+x[1]**2

x,varnames=dakota_utils.read_input(sys.argv[1])[0]

fval = my_func(x)

dakota_utils.write_output(sys.argv[2],fval,'my_function')



    
