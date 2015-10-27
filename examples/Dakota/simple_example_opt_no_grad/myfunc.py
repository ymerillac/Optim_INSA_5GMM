import sys
import numpy
#import dakota_utils #--> le module a creer

def my_func(x):
    return x[0]**2+x[1]**2

inp_file = sys.argv[1]

fid = open(inp_file,'r')

#print fid.read() #--> pour voir les details

# read variables
# x,i = dakota_utils.read_input(sys.argv[1]) #--> notre but final
first_line = fid.readline()
n_var = eval(first_line.split()[0])
#print "number of variables : ",n_var #--> afficher le nombre de variables
# read variables values
x = numpy.zeros(n_var)
varnames = []
for i in xrange(n_var):
    line = fid.readline()
    words = line.split()
    x[i] = eval(words[0])
    varnames.append(words[1].strip())
while "eval_id" not in line:
    line=fid.readline()
line_tmp=line.split(" ")
nit=line_tmp[-2]
print "iteration : ",nit
fid.close()

fval = my_func(x)
print "valeur de la fonction :",fval
print ("arguments : "+varnames[0]+" = {}, ".format(x[0])+varnames[1]+" = {}\n".format(x[1]))

# write output file
fid = open(sys.argv[2],'w')
fid.write(str(fval)+'    my_function\n')
fid.close()




    
