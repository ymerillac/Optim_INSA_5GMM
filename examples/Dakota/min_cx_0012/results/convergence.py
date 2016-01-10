import sys
from matplotlib.pyplot import *
import numpy

fichier_genetique=open('optim_res_genetique.dat','r')
data=fichier_genetique.readlines()
tmp=data[-1]
tmp=tmp.split()
nbit_genetique=int(tmp[0])
cost_function_genetique=numpy.zeros(nbit_genetique)
data=data[1:nbit_genetique+1]
nbok=0
for i in xrange(nbit_genetique):
	tmp=data[i]
	tmp=tmp.split()
	cstr=float(tmp[-1])
	if cstr<=0:
		cost_function_genetique[nbok]=float(tmp[-2]) # -2 car my_constraint
		nbok += 1
cost_function_genetique=cost_function_genetique[:nbok]
fichier_genetique.close()

fichier_quasi_newton=open('optim_res_quasi_newton.dat','r')
data=fichier_quasi_newton.readlines()
tmp=data[-1]
tmp=tmp.split()
nbit_quasi_newton=int(tmp[0])
cost_function_quasi_newton=numpy.zeros(nbit_quasi_newton)
data=data[1:nbit_quasi_newton+1]
nbok2=0
for i in xrange(nbit_quasi_newton):
	tmp=data[i]
	tmp=tmp.split()
	cstr=float(tmp[-1])
	if cstr<=0:
		cost_function_quasi_newton[nbok2]=float(tmp[-2]) # -2 car my_constraint
		nbok2 += 1
cost_function_quasi_newton=cost_function_quasi_newton[:nbok2]
fichier_quasi_newton.close()

print 'Quasi newton :',nbit_quasi_newton,'iterations. Genetique :',nbit_genetique,'iterations.'
print 'Final cost function value'
print 'Quasi newton :',cost_function_quasi_newton[-1],' Genetique : ',cost_function_genetique[-1]

figure()
plot(range(1,nbok+1),cost_function_genetique,'+',label='genetique')
title("fonction cout")
plot(cost_function_quasi_newton,'r',label='quasi newton')
legend()
show()


