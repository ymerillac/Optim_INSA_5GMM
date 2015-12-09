import sys
import numpy
import dakota_utils #--> le module a creer

def my_func(x):
    return x[0]**2+x[1]**2

def constraint(x):
	return x[0]-2.

x,varnames=dakota_utils.read_input(sys.argv[1])

fval = my_func(x)
cstr = constraint(x)

dakota_utils.write_outputs(sys.argv[2],[fval,cstr],['my_function','my_constraint'])



    
