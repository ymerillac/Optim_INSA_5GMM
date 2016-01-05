import NACA 
import numpy 
import scipy.optimize as scop 
import build_mesh
from apame_utils import run_case
import mesh_utils
from scipy import interpolate

nb_iter = 500
cx_tab = numpy.array(500)
contrainte_tab = numpy.array(500)
k=0

def my_func(X):
	chord=[1.,1.]
	M=X[0:2]
	P=X[2:4]
	T=X[4:6]
	alpha=X[6]
	span=X[7]
	chord[1]=X[8]
	n_sections=31
	n_naca_pts=50
	Sref=10.
	Re=5.e5
	mesh_utils.write_wing_file(M,P,T,chord,span,n_sections,n_naca_pts)
	NACA.create_wing('current_wing','output')
	run_apame(alpha)
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
	#f = interpolate.interp1d(numpy.array(cz), numpy.array(cx), kind='cubic')
	#print 'DRAG >>> ', min(cx)
	#cx=f(0.5)
	
	ind = cx.argmin()
	Cf = 0.0583/Re**0.2
	Sw = (chord[0]+chord[1])*span/2.
	Cd_friction = Sw/Sref*Cf
	print 'DRAG >>> ', cx+Cd_friction
	cx_tab[k]=cx+Cd_friction
	contrainte_tab[k]=0.5-cz
	k=k+1
	return cx+Cd_friction#,cz[ind]

def run_apame(alpha):
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
	#nrun=11
	#alpha_list = numpy.linspace(-10.,10.,nrun)
	alpha_list = numpy.array([alpha])
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
		 Sref=10.,
		 method=0,
		 farfield_dist=50.,
		 velorder=1)


## Test de scipy 
x0 = numpy.array([0.]*2+[0.4]*2+[0.12]*2+[5.78]*1+[10.]*1+[1.]*1)
bounds = [(-0.05,0.05)]*2+[(0.3,0.5)]*2+[(0.08,0.15)]*2+[(-10.,10.)]*1+[(5.,15.)]*1+[(0.5,1.5)]*1
[res,v=scop.fmin_l_bfgs_b(my_func,x0,approx_grad=True,bounds=bounds,maxfun=nb_iter)
#my_func(x0)


