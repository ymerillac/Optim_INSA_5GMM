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

def create_mesh_linear_interp(m_root,p_root,t_root,m_tip,p_tip,t_tip,ystart,yend,chord,n_sections,n_naca_points=150):
    if n_sections%2 == 0:
	raise ValueError("Merci de donner un nombre de sections impair")
    ind_root=(n_sections-1)//2
    semi_span = (ystart-yend)/2
    theta = numpy.linspace(0.,numpy.pi,n_sections)
    y_list = -semi_span*numpy.cos(theta)
    y_list[ind_root]=0. # souvent 10^-16

    Xu_tip,Xl_tip,Yu_tip,Yl_tip = create(m_tip,p_tip,t_tip,chord,n_naca_points)
    upper_tip = numpy.array([Xu_tip,Yu_tip]).T
    lower_tip = numpy.array([Xl_tip,Yl_tip]).T
    tip_section = numpy.concatenate((upper_tip[::-1],lower_tip[1:]))

    Xu_root,Xl_root,Yu_root,Yl_root = create(m_root,p_root,t_root,chord,n_naca_points)
    upper_root = numpy.array([Xu_root,Yu_root]).T
    lower_root = numpy.array([Xl_root,Yl_root]).T
    root_section = numpy.concatenate((upper_root[::-1],lower_root[1:]))

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
    writer.SetFileName('output.vtp')
    writer.Write()


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

create_mesh_linear_interp(0.,0.12,0.4,0.1,0.15,0.3,-2.5,2.5,1.,21,50)
