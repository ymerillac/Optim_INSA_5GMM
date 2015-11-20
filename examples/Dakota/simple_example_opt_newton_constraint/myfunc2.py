import sys
import numpy
import math
import dakota_utils #--> le module a creer

def my_func(x):
    return math.exp(x[0]**2+x[1]**2)

x,varnames=dakota_utils.read_input(sys.argv[1])

fval = my_func(x)

dakota_utils.write_output(sys.argv[2],fval)



