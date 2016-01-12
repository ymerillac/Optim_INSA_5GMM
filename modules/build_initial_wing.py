import sys
import mesh_utils
import NACA
from apame_utils import run_case
import numpy

nb=int(sys.argv[1])
n_sections=31
n_naca_pts=50
chord=[1.,1.]
span=10.
Re=5.e5
Sref=10.

def run_apame(nb_angles):
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
	alpha_span=float((nb_angles-1)//2)
	alpha_list = numpy.linspace(-alpha_span,alpha_span,nb_angles)
	beta_list = numpy.zeros(nb_angles)

	run_case('initial'+str(nb)+'.vtp',
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

if nb==8415:
	nb_angles=11
	M=[0.08,0.08]
	P=[0.4,0.4]
	T=[0.15,0.15]
	
elif nb==12:
	nb_angles=21
	M=[0.0,0.0]
	P=[0.4,0.4]
	T=[0.12,0.12]

mesh_utils.write_wing_file(M,P,T,chord,span,n_sections,n_naca_pts)
NACA.create_wing('current_wing','initial'+str(nb))
run_apame(nb_angles)


