from NACA import create_wing
import numpy
import vtk
import dakota_utils

def write_test_file(m,p,t,ystart,yend,nsections,nnaca):
    #y_list = numpy.linspace(ystart,yend,nsections)
    theta=numpy.linspace(0.,numpy.pi,nsections)
    span=yend-ystart
    y_list=-0.5*span*numpy.cos(theta)
    fid = open('wing','w')
    fid.write("# id variable              valeur\n")
    fid.write("PROFILE.N_SECTIONS "+str(nsections)+"\n")
    fid.write("PROFILE.N_NACA_PTS "+str(nnaca)+"\n")
    i=0
    fid.write('SECTION_ROOT.Y '+str(0.)+'\n')
    fid.write('SECTION_ROOT.CORDE '+'1.0\n')
    fid.write('SECTION_ROOT.M '+str(m[i])+'\n')
    fid.write('SECTION_ROOT.P '+str(p[i])+'\n')
    fid.write('SECTION_ROOT.T '+str(t[i])+'\n')
    i=1
    fid.write('SECTION_TIP.Y '+str(-5.)+'\n')
    fid.write('SECTION_TIP.CORDE '+'1.0\n')
    fid.write('SECTION_TIP.M '+str(m[i])+'\n')
    fid.write('SECTION_TIP.P '+str(p[i])+'\n')
    fid.write('SECTION_TIP.T '+str(t[i])+'\n')

    fid.close()

def create_mesh(M,P,T,nsections):
	n_naca_pts = 50
	write_test_file(M,P,T,-5.,5.,nsections,n_naca_pts)
	create_wing("wing","output")

