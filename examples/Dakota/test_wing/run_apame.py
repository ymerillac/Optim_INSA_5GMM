from apame_utils import run_case
import numpy

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
alpha_list = numpy.linspace(-10.,10.,10)
beta_list = numpy.zeros(10)

run_case('naca0012.vtp',
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
