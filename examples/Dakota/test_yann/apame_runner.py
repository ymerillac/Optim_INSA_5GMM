import sys
import numpy
import dakota_utils
from build_mesh import build_mesh
from apame_utils import run_case

x,varnames,=dakota_utils.read_input(sys.argv[1])

filename = 'wing.vtp'

build_mesh(filename,
           semi_span=x[0],
           root_chord=x[1],
           tip_chord=x[2],
           root_m=x[3],
           root_p=x[4],
           tip_m=x[5],
           tip_p=x[6],
           t=0.1,
           nu=50,
           nv=21)

Patm = 101325.#Pa
rho = 1.225#kg.m-3
c = numpy.sqrt(1.4*Patm/rho)
Mach = 0.17
airspeed = c*Mach

run_case(filename,
         wake_length=20.,
         alpha=[x[7]],
         beta=[0.],
         v=25.,
         rho=rho,
         P=Patm,
         Mach=Mach,
         origin=[0.,0.,0.],
         wingspan=2.*x[0],
         ref_chord=1.,
         Sref=10.,
         method=0,
         farfield_dist=50.,
         velorder=1)

fid = open('wing_polar.dat')
lines = fid.readlines()
fid.close()
[_,Cl,Cd] = [eval(word) for word in lines[0].split()]
    
print Cl,10000.*Cd

dakota_utils.write_outputs(sys.argv[2],[Cd,0.5-Cl],['Cd','Cl'])