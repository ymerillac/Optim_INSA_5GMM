from NACA import create_wing
from NACA import create
import numpy
import vtk
import dakota_utils
import matplotlib.pyplot as plt

def write_test_file(m,p,t,ystart,yend,nsections=5):
    y_list = numpy.linspace(ystart,yend,nsections)
    fid = open('wing','w')
    fid.write("# id variable              valeur\n")
    for i in xrange(nsections):
        fid.write('SECTION'+str(i+1)+'.Y     '+str(y_list[i])+'\n')
        fid.write('SECTION'+str(i+1)+'.CORDE 1.0\n')
        fid.write('SECTION'+str(i+1)+'.M     '+str(m[i])+'\n')
        fid.write('SECTION'+str(i+1)+'.P     '+str(p[i])+'\n')
        fid.write('SECTION'+str(i+1)+'.T     '+str(t[i])+'\n')
    fid.close()

def create_profile_linear_interp(m_root,p_root,t_root,m_tip,p_tip,t_tip,ystart,yend,chord,n_sections,n_naca_points=150):
    if n_sections%2 == 0:
	raise ValueError("Merci de donner un nombre de sections impair")
    ind_root=(n_sections-1)//2
    semi_span = (ystart-yend)/2
    theta = numpy.linspace(0.,numpy.pi,n_sections)
    y_list = -semi_span*numpy.cos(theta)
    y_list[ind_root]=0. # souvent 10^-16

    tip_section_begin = list(create(m_root,p_root,t_root,chord,n_naca_points))
    root_section = list(create(m_tip,p_tip,t_tip,chord,n_naca_points))
    tip_section_end = list(create(m_root,p_root,t_root,chord,n_naca_points))

    tip_section_begin.insert(0,-semi_span)
    root_section.insert(0,y_list[ind_root])
    tip_section_end.insert(0,semi_span)

    profile = [tip_section_begin]

    for i in xrange(1,n_sections-1):
        if i == ind_root:
	    profile.append(root_section)	    
	else:
	    r = numpy.abs(y_list[i])/semi_span
	    current_section = [y_list[i],(1.-r)*root_section[1] + r*tip_section_begin[1], \
            (1.-r)*root_section[2] + r*tip_section_begin[2],(1.-r)*root_section[3] + r*tip_section_begin[3], \
            (1.-r)*root_section[4] + r*tip_section_begin[4]]
	    profile.append(current_section)
    profile.append(tip_section_end)

    return profile

    #ch = str(chord)
    #fid = open('wing','w')
    #fid.write("# id variable              valeur\n")
    #for i in xrange(n_sections):
     #   fid.write('SECTION'+str(i+1)+'.Y     '+str(y_list[i])+'\n')
      #  fid.write('SECTION'+str(i+1)+'.CORDE '+ch+'\n')
       # fid.write('SECTION'+str(i+1)+'.M     '+str(m[i])+'\n')
        #fid.write('SECTION'+str(i+1)+'.P     '+str(p[i])+'\n')
        #fid.write('SECTION'+str(i+1)+'.T     '+str(t[i])+'\n')
   # fid.close()

    return profile

def create_mesh(M,P,T):
	n_naca_pts = 50
	write_test_file(M,P,T,-5.,5.,nsections=5)
	wing=create_wing("wing",n_naca_pts)
	n_sections=len(wing)
	n_section_pts = 2*n_naca_pts-1

	# build mesh
	vtk_model = vtk.vtkStructuredGrid()
	vtk_model.SetDimensions(n_section_pts,n_sections,1)
	# build points
	vtk_points = vtk.vtkPoints()

	for j in xrange(n_sections):
	    upper_pts = numpy.array([wing[j][1],wing[j][3]]).T
	    lower_pts = numpy.array([wing[j][2],wing[j][4]]).T
	    section_pts = numpy.concatenate((lower_pts[::-1],upper_pts[1:]))
	    for i in xrange(n_section_pts):
		vtk_points.InsertNextPoint(section_pts[i,0],wing[j][0],section_pts[i,1])
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
	writer.SetFileName('output.vtp')
	writer.Write()

def create_mesh_linear_interp(m_root,p_root,t_root,m_tip,p_tip,t_tip,ystart,yend,chord,n_sections):
	n_naca_pts = 50
	wing=create_profile_linear_interp(m_root,p_root,t_root,m_tip,p_tip,t_tip,ystart,yend,chord,n_sections,n_naca_pts)
	n_section_pts = 2*n_naca_pts-1

	# build mesh
	vtk_model = vtk.vtkStructuredGrid()
	vtk_model.SetDimensions(n_section_pts,n_sections,1)
	# build points
	vtk_points = vtk.vtkPoints()

	for j in xrange(n_sections):
	    upper_pts = numpy.array([wing[j][1],wing[j][3]]).T
	    lower_pts = numpy.array([wing[j][2],wing[j][4]]).T
	    section_pts = numpy.concatenate((lower_pts[::-1],upper_pts[1:]))
	    for i in xrange(n_section_pts):
		vtk_points.InsertNextPoint(section_pts[i,0],wing[j][0],section_pts[i,1])
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
	writer.SetFileName('output.vtp')
	writer.Write()

create_mesh_linear_interp(0.,0.4,0.12,0.01,0.38,0.10,-5.,5.,1.,15)
