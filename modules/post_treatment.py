import sys
from matplotlib.pyplot import *
import numpy
from math import sqrt

inp_file=sys.argv[1] 
fichier=open(inp_file,'r')
nb_it=0
first_line=fichier.readline()
first_line=first_line.split()
first_line=first_line[2:-1]
nvar=len(first_line)
while fichier.readline():
	nb_it+=1
fichier.close()
var=numpy.zeros((nb_it,nvar))
func=numpy.zeros(nb_it)	
fichier=open(inp_file,'r')
fichier.readline()
for l in xrange(nb_it):
	tmp=fichier.readline()
	tmp=tmp.split()
	func[l]=float(tmp[-1])
	tmp=[float(s) for s in tmp[2:-1]]
	var[l][:]=tmp
fichier.close()
figure()
plot(range(1,nb_it+1),func)
title("fonction cout")
sq=sqrt(nvar)
for col in xrange(nvar):
	figure()
	plot(range(1,nb_it+1),var.transpose()[:][col])
	title(first_line[col])
show()


