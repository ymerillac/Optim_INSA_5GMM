import sys
import numpy

def my_func(x):
    return x[0]**2+x[1]**2

inp_file = sys.argv[1]

fid = open(inp_file,'r')

# read variables
first_line = fid.readline()
n_var = eval(first_line.split()[0])
print "number of variables : ",n_var
# read variables values
x = numpy.zeros(n_var)
varnames = []
for i in xrange(n_var):
    line = fid.readline()
    words = line.split()
    x[i] = eval(words[0])
    varnames.append(words[1].strip())
fid.close()

fval = my_func(x)

# write output file
fid = open(sys.argv[2],'w')
fid.write(str(fval)+'    my_function\n')
fid.close()




    
