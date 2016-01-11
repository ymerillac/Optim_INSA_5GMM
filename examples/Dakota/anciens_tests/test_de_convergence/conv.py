import numpy
from NACA import create_wing
import vtk
import dakota_utils
import matplotlib.pyplot as plt
from build_mesh import create_mesh
from run_apame import run

N_sections=numpy.linspace(5,101,49)
Cx=numpy.zeros(49)
Cz=numpy.zeros(49)
it=0
for i in N_sections:
	print('Test de convergence - section numero {}'.format(int(i)))
	M=0.08*numpy.ones(2)
	P=0.4*numpy.ones(2)
	T=0.15*numpy.ones(2)
	create_mesh(M,P,T,int(i))
	run()
	fid = open('output_polar.dat','r')
	tmp=fid.readline()
	tmp=tmp.split()
	Cx[it]=float(tmp[2])
	Cz[it]=float(tmp[1])
	it+=1
	fid.close()

plt.figure()
plt.plot(N_sections,Cx)
plt.title('Cx')
plt.figure()
plt.plot(N_sections,Cz)
plt.title('Cz')
plt.show()
