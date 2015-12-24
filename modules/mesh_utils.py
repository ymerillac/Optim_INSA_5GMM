import NACA
import numpy
import vtk
import dakota_utils
import matplotlib.pyplot as plt



def create_mesh_linear_interp(filename,root,tip,span,n_sections,n_naca_points=150):
    if n_sections%2 == 0:
	raise ValueError("Merci de donner un nombre de sections impair")
    if type(root) is not list:
	raise TypeError("Le deuxieme argument doit etre une liste contenant M,P,T et Chord pour la section root.")
    if type(tip) is not list:
	raise TypeError("Le troisieme argument doit etre une liste contenant M,P,T et Chord pour la section tip.")
    if type(filename) is not str:
	raise TypeError("Le premier argument doit etre un string (nom du fichier de sortie).")
    if len(root)!=4:
	raise ValueError("Le deuxieme argument doit etre une liste contenant 4 valeurs (M,P,T,Chord).")
    if len(tip)!=4:
	raise ValueError("Le troisieme argument doit etre une liste contenant 4 valeurs (M,P,T,Chord).")

    ind_root=(n_sections-1)//2
    semi_span = span/2.
    theta = numpy.linspace(0.,numpy.pi,n_sections)
    y_list = -semi_span*numpy.cos(theta)
    y_list[ind_root]=0. # souvent 10^-16
    m_root,p_root,t_root,chord_root=root
    m_tip,p_tip,t_tip,chord_tip=tip

    Xu_tip,Xl_tip,Yu_tip,Yl_tip = NACA.create(m_tip,p_tip,t_tip,chord_tip,n_naca_points)
    upper_tip = numpy.array([Xu_tip,Yu_tip]).T
    lower_tip = numpy.array([Xl_tip,Yl_tip]).T
    tip_section = numpy.concatenate((lower_tip[::-1],upper_tip[1:]))

    Xu_root,Xl_root,Yu_root,Yl_root = NACA.create(m_root,p_root,t_root,chord_root,n_naca_points)
    upper_root = numpy.array([Xu_root,Yu_root]).T
    lower_root = numpy.array([Xl_root,Yl_root]).T
    root_section = numpy.concatenate((lower_root[::-1],upper_root[1:]))


    # build mesh
    vtk_model = vtk.vtkStructuredGrid()
    vtk_model.SetDimensions(2*n_naca_points-1,n_sections,1)
    # build points
    vtk_points = vtk.vtkPoints()

    for i in xrange(n_sections):
	r = numpy.abs(y_list[i])/semi_span
	current_section = (1.-r)*root_section + r*tip_section
        for j in xrange(2*n_naca_points-1):
                vtk_points.InsertNextPoint(current_section[j,0],y_list[i],current_section[j,1])
    # set points
    vtk_model.SetPoints(vtk_points)
    # convert to poly data    
    pdata_filter = vtk.vtkGeometryFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        pdata_filter.SetInput(vtk_model)
    else:
        pdata_filter.SetInputData(vtk_model)
    pdata_filter.Update()
    poly_data = pdata_filter.GetOutput()
    
    # compute normals
    norms = vtk.vtkPolyDataNormals()
    if vtk.VTK_MAJOR_VERSION <= 5:
        norms.SetInput(poly_data)
    else:
        norms.SetInputData(poly_data)
    norms.ComputePointNormalsOff()
    norms.ComputeCellNormalsOn()
    norms.ConsistencyOn()
    norms.Update()
    
    # clean poly data
    clean_poly = vtk.vtkCleanPolyData()
    clean_poly.ToleranceIsAbsoluteOn()
    clean_poly.SetAbsoluteTolerance(1.e-6)
    if vtk.VTK_MAJOR_VERSION <= 5:
        clean_poly.SetInput(norms.GetOutput())
    else:
        clean_poly.SetInputData(norms.GetOutput())
    clean_poly.Update()
    
    # write output mesh
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInput(clean_poly.GetOutput())
    else:
        writer.SetInputData(clean_poly.GetOutput())
    writer.SetFileName(filename+'.vtp')
    writer.Write()



def write_wing_file(m,p,t,chord,span,n_sections,n_naca_points):
    # m,p,t,chord des listes contenant respectivement les valeurs de m_root et m_tip, p_root et p_tip etc.
    # a donner a manger a myfunc
    fid = open('current_wing','w')
    fid.write("variable                valeur\n")
    fid.write('SECTION_ROOT.Y         '+str(0.)+'\n')
    fid.write('SECTION_ROOT.CORDE     '+str(chord[0])+'\n')
    fid.write('SECTION_ROOT.M         '+str(m[0])+'\n')
    fid.write('SECTION_ROOT.P         '+str(p[0])+'\n')
    fid.write('SECTION_ROOT.T         '+str(t[0])+'\n')
    fid.write('SECTION_TIP.Y          '+str(-span/2.)+'\n')
    fid.write('SECTION_TIP.CORDE      '+str(chord[1])+'\n')
    fid.write('SECTION_TIP.M          '+str(m[1])+'\n')
    fid.write('SECTION_TIP.P          '+str(p[1])+'\n')
    fid.write('SECTION_TIP.T          '+str(t[1])+'\n')
    fid.write('PROFILE.N_SECTIONS     '+str(n_sections)+'\n')
    fid.write('PROFILE.N_NACA_PTS     '+str(n_naca_points)+'\n')
    fid.close()
