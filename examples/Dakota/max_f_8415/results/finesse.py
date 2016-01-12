import sys
from matplotlib.pyplot import *
import numpy

# ordre d'apparition des argv : initial, bfgs, quasi-newton

Sref = 10.
Re = 5.e5
Cf = 0.0583/Re**0.2
chord = [1.,1.]
alpha = numpy.linspace(-5.,5.,11)
cz = numpy.zeros(11)
cx = numpy.zeros(11)
finesse_initial = numpy.zeros(11)
finesse_bfgs = numpy.zeros(11)
finesse_qn = numpy.zeros(11)
liste_parametres = [sys.argv[1],sys.argv[3],sys.argv[5]]
liste_donnees = [sys.argv[2],sys.argv[4],sys.argv[6]]
liste_finesse = [finesse_initial,finesse_bfgs,finesse_qn]

for j in xrange(3):
	fichier_parametres = open(liste_parametres[j],'r')
	lines = fichier_parametres.readlines()
	fichier_parametres.close()
	lines = lines[-1]
	lines = lines.split()
	chord[1] = float(lines[9])
	span = float(lines[8])
	Sw = (chord[0]+chord[1])*span/2.
	Cd_friction = Sw/Sref*Cf

	fichier_donnees = open(liste_donnees[j],'r')
	donnees = fichier_donnees.readlines()
	fichier_donnees.close()
	for i,tmp in enumerate(donnees):
		tmp = tmp.split()
		print tmp
		cz[i] = float(tmp[1])
		cx[i] = float(tmp[2])+Cd_friction
	liste_finesse[j] = cz/cx

figure()
plot(alpha,liste_finesse[0],color="k",label="Initial")
plot(alpha,liste_finesse[1],color="b",label="BFGS")
plot(alpha,liste_finesse[2],color="r",label="Quasi-Newton")
title("L/D ratio")
legend(loc=4)
xlabel("alpha")
show()

	
