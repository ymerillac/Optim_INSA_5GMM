from apame_utils import run_case
import math
import numpy

gamma = 1.4
Patm = 101325.#Pa
rho = 1.225#kg.m-3
Mach = 0.5
c = math.sqrt(gamma*Patm/rho)
mu = 1.5e-5
Airspeed = c*Mach
origin=[0.,0.,0.]
wingspan=2.3926
ref_chord=0.8059
Sref=1.
Re = rho*Airspeed*ref_chord/mu
Cf = 0.0583/Re**0.2
Sw = (ref_chord)*wingspan
Cd_friction = Sw/Sref*Cf
#run_case('onera_m6wing.vtp',10.,numpy.linspace(-5.,5.,11),numpy.zeros(11),Airspeed,rho,Patm,Mach,origin,wingspan,ref_chord,Sref,0,5.,1)

fichier_donnees = open('polar_m6.txt')
donnees_cx = numpy.zeros(11)
donnees_cz = numpy.zeros(11)
tmp = fichier_donnees.readlines()
tmp = tmp[1:]
for i in xrange(11):
	tmp2 = tmp[i].split()
	tmp_cz = tmp2[-2]
	tmp_cx = tmp2[-1]
	donnees_cx[i] = float(tmp_cx)
	donnees_cz[i] = float(tmp_cz)
fichier_donnees.close()

fichier_apame = open('onera_m6wing_polar.dat')
apame_cx = numpy.zeros(11)
apame_cz = numpy.zeros(11)
tmp = fichier_apame.readlines()
for i in xrange(11):
	tmp2 = tmp[i].split()
	tmp_cz = tmp2[-2]
	tmp_cx = tmp2[-1]
	apame_cx[i] = float(tmp_cx) #+ Cd_friction
	apame_cz[i] = float(tmp_cz)
fichier_apame.close()

diff_cx = numpy.zeros(11)
diff_cx = 100.*numpy.divide(abs(apame_cx-donnees_cx),abs(donnees_cx))
diff_cz = numpy.zeros(11)
diff_cz = 100.*numpy.divide(abs(apame_cz-donnees_cz),abs(donnees_cz))
print 'Ecart relatif, Cx\n',diff_cx
print '\nEcart relatif, Cz\n',diff_cz

