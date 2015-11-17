from NACA import create as create_naca
import numpy
import vtk

semi_span = 5.
n_sections = 10
n_naca_pts = 50
n_section_pts = 2*n_naca_pts-1
chord = 1.
m = 0.
p = 0.4
t = 0.12
y_sections = numpy.linspace(-semi_span,semi_span,n_sections)

# build mesh
vtk_model = vtk.vtkStructuredGrid()
vtk_model.SetDimensions(n_section_pts,n_sections,1)
# build points
vtk_points = vtk.vtkPoints()
for i,y in enumerate(y_sections):
    Xu,Xl,Yu,Yl=create_naca(m,p,t,C=chord,n=n_naca_pts)
    upper_pts = numpy.array([Xu,Yu]).T
    lower_pts = numpy.array([Xl,Yl]).T
    section_pts = numpy.concatenate((lower_pts[::-1],upper_pts[1:]))
    for j in xrange(n_section_pts):
        vtk_points.InsertNextPoint(section_pts[j,0],y,section_pts[j,1])
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
writer.SetFileName('naca0012.vtp')
writer.Write()


