import sys
from matplotlib.pyplot import *
import numpy
from math import sqrt
from math import ceil

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
for l in xrange(nb_it-1): # on suppose que dans la derniere ligne on a deja recopie la meilleure valeur
	tmp=fichier.readline()
	tmp=tmp.split()
	func[l]=float(tmp[-1])
	tmp=[float(s) for s in tmp[2:-1]]
	var[l][:]=tmp
fichier.close()
figure()
plot(range(1,nb_it+1),func)
title("fonction cout")
nb_subplot=nvar//4
last_plot=nb_subplot%4
for nb_fig in xrange(nb_subplot):
	figure()
	subplot(221)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_fig)*4])
	title(first_line[(nb_fig)*4])
	subplot(222)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_fig)*4+1])
	title(first_line[(nb_fig)*4+1])
	subplot(223)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_fig)*4+2])
	title(first_line[(nb_fig)*4+2])
	subplot(224)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_fig)*4+3])
	title(first_line[(nb_fig)*4+3])
if last_plot==3:
	figure()
	subplot(221)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_subplot)*4])
	title(first_line[(nb_fig)*4])
	subplot(222)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_subplot)*4+1])
	title(first_line[(nb_fig)*4+1])
	subplot(223)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_subplot)*4+2])
	title(first_line[(nb_fig)*4+2])
elif last_plot==2:
	figure()
	subplot(121)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_subplot)*4])
	title(first_line[(nb_fig)*4])
	subplot(122)
	plot(range(1,nb_it+1),var.transpose()[:][(nb_subplot)*4+1])
	title(first_line[(nb_fig)*4+1])
elif last_plot==1:
	figure()
	plot(range(1,nb_it+1),var.transpose()[:][(nb_subplot)*4])
	title(first_line[(nb_fig)*4])
show()


