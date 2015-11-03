import sys
import numpy

def read_input(inp_file):
    fichier=open(inp_file,'r')
    first_line=fichier.readline()
    n_var = eval(first_line.split()[0])
    x = numpy.zeros(n_var)
    varnames = []
    for i in xrange(n_var):
    	line = fichier.readline()
        words = line.split()
        x[i] = eval(words[0])
        varnames.append(words[1].strip())
    fichier.close()
    return x

def write_output(out_file,fval,name):
    fichier=open(sys.argv[2],'w')
    fichier.write(str(fval)+'    '+str(name)+'\n')
    fichier.close()

