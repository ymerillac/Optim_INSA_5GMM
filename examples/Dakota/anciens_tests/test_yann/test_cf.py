import sys
import numpy
from build_mesh import build_mesh
from apame_utils import run_case

semi_span=5.*1.1
root_chord=2./1.1
tip_chord=1./1.1
root_m=0.01
root_p=0.4
tip_m=0.01
tip_p=0.4
t = 0.11


filename = 'wing.vtp'

build_mesh(filename,
           semi_span,
           root_chord,
           tip_chord,
           root_m,
           root_p,
           tip_m,
           tip_p,
           t=t,
           nu=50,
           nv=41)

Patm = 101325.#Pa
rho = 1.225#kg.m-3
c = numpy.sqrt(1.4*Patm/rho)
Mach = 0.15
airspeed = c*Mach
mu=1.5e-5
Re = rho*airspeed*root_chord/mu
print "Re=",Re
print "airspeed=",airspeed*3.6

Sref=15.
Sw = (root_chord+tip_chord)*semi_span

print "Sw=",Sw

run_case(filename,
         wake_length=20.,
         alpha=[2.],
         beta=[0.],
         v=airspeed,
         rho=rho,
         P=Patm,
         Mach=Mach,
         origin=[0.,0.,0.],
         wingspan=2.*semi_span,
         ref_chord=root_chord,
         Sref=Sref,
         method=0,
         farfield_dist=50.,
         velorder=1)

fid = open('wing_polar.dat')
lines = fid.readlines()
fid.close()
[_,Cl,Cd] = [eval(word) for word in lines[0].split()]
    
print Cl,10000.*Cd

Qd = 0.5*rho*Sref*airspeed**2

Cf=0.0583/Re**0.2

print "Lift force (N) : ",Qd*Cl
print "Pressure Drag force (N) : ",Qd*Cd
print "Friction Drag force (N) : ",0.5*rho*Sw*airspeed**2*Cf

print '**** drag coefficients *****'
print "Pressure: ",Cd
print "Friction: ",Sw/Sref*Cf


