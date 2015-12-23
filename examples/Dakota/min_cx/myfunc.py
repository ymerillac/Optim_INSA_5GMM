import sys
import numpy
import dakota_utils
import mesh_utils
import NACA
from apame_utils import run_case

# constantes
n_sections=31
n_naca_pts=50
chord=[1.,1.]
Sref=10.
Re=5.e5

def my_func(M,P,T,aoa,span,tip_chord):
        chord[1] = tip_chord
	mesh_utils.write_wing_file(M,P,T,chord,span,n_sections,n_naca_pts)
	NACA.create_wing('current_wing','output')
	run_apame(aoa)
	fid = open('output_polar.dat','r')
	nb_lines=0
	while fid.readline():
		nb_lines+=1
	if nb_lines==0:
		raise TypeError("Empty file")
	fid.close()
	fid = open('output_polar.dat','r')
	# revenir au debut
	alpha=numpy.zeros(nb_lines)
	cx=numpy.zeros(nb_lines)
	cz=numpy.zeros(nb_lines)
	for i in xrange(nb_lines):
		tmp=fid.readline()
		tmp=tmp.split()
		alpha[i]=float(tmp[0])
		cx[i]=float(tmp[2])
		cz[i]=float(tmp[1])
	fid.close()
	ind = cx.argmin()
	Cf = 0.0583/Re**0.2
	Sw = (chord[0]+chord[1])*span/2.
        Cd_friction = Sw/Sref*Cf
	return min(cx)+Cd_friction,cz[ind]


def run_apame(aoa):
	# compute flight conditions
	r=287.058 #J.kg-1.K-1
	Patm = 101325.#Pa
	rho = 1.225#kg.m-3
	mu=1.5e-5
	T = Patm/(r*rho)
	Tc = T-273.15
	print "T=",T

	airspeed = Re*mu/rho

	print "Airspeed=",airspeed*3.6

	# list of alpha and beta to compute
	alpha_list = numpy.array([aoa])
	beta_list = numpy.zeros(1)

	run_case('output.vtp',
		 wake_length=10.,
		 alpha=alpha_list,
		 beta=beta_list,
		 v=airspeed,
		 rho=rho,
		 P=Patm,
		 Mach=0.,
		 origin=[0.,0.,0.],
		 wingspan=10.,
		 ref_chord=1.,
		 Sref=Sref,
		 method=0,
		 farfield_dist=50.,
		 velorder=1)

X,varnames=dakota_utils.read_input(sys.argv[1])
M = [X[0],X[1]]
P = [X[2],X[3]]
T = [X[4],X[5]]
aoa = X[6]
span = X[7]
tip_chord = X[8]

fval,cstr = my_func(M,P,T,aoa,span,tip_chord)
dakota_utils.write_outputs(sys.argv[2],[fval,-cstr+0.5],['my_function','my_constraint'])
    
