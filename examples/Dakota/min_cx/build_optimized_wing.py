import mesh_utils
import numpy
from apame_utils import run_case

def run_apame():
	# compute flight conditions
	r=287.058 #J.kg-1.K-1
	Patm = 101325.#Pa
	rho = 1.225#kg.m-3
	mu=1.5e-5
	T = Patm/(r*rho)
	Tc = T-273.15
	print "T=",T

	Re=5.e5
	airspeed = Re*mu/rho

	print "Airspeed=",airspeed*3.6

	# list of alpha and beta to compute
	alpha_list = numpy.array([3])
	beta_list = numpy.zeros(1)

	run_case('optimized_wing.vtp',
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
		 Sref=10.,
		 method=0,
		 farfield_dist=50.,
		 velorder=1)

file = open('optim_res.dat','r')
tmp = file.readlines()
tmp = tmp[-1]
tmp = tmp.split()
print tmp
# attention : a changer quand il y aura plus de variables de design. Dernier argument : la corde
root = [float(tmp[2]),float(tmp[4]),float(tmp[6]),1.]
tip = [float(tmp[3]),float(tmp[5]),float(tmp[7]),1.]
mesh_utils.create_mesh_linear_interp('optimized_wing',root,tip,10.,31,50)
run_apame()
file.close()
