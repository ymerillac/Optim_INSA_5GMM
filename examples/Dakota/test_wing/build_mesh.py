from NACA import create_wing
import numpy
import vtk
import dakota_utils

def write_test_file(m,p,t,ystart,yend,nsections=5):
    y_list = numpy.linspace(ystart,yend,nsections)
    fid = open('TEST','w')
    fid.write("# id variable              valeur\n")
    for i in xrange(nsections):
        fid.write('SECTION'+str(i+1)+'.Y     '+str(y_list[i])+'\n')
        fid.write('SECTION'+str(i+1)+'.CORDE 1.0\n')
        fid.write('SECTION'+str(i+1)+'.M     '+str(m[i])+'\n')
        fid.write('SECTION'+str(i+1)+'.P     '+str(p[i])+'\n')
        fid.write('SECTION'+str(i+1)+'.T     '+str(t[i])+'\n')
    fid.close()

n_naca_pts = 50
M=[0.,0.1,0.1,0.4,0.7]
P=[0.4, 0.4, 0.4, 0.4, 0.4] 
T=[0.12,0.12,0.12,0.12,0.14] 
write_test_file(M,P,T,-5.,5.,nsections=5)
wing=create_wing("TEST",n_naca_pts)
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
writer.SetFileName('naca0012.vtp')
writer.Write()


