from NACA import create as create_naca
import numpy
import vtk


def build_mesh(filename,semi_span,root_chord,tip_chord,root_m,root_p,tip_m,tip_p,t,nu=50,nv=21):
    # build root section
    Xu,Xl,Yu,Yl = create_naca(root_m,root_p,t,C=root_chord,n=nu)
    upper = numpy.array([Xu,Yu]).T
    lower = numpy.array([Xl,Yl]).T
    root_section = numpy.concatenate((lower[::-1],upper[1:]))
    # set x position
    root_section[:,0]-=0.25*root_chord
    # build tip section
    Xu,Xl,Yu,Yl = create_naca(tip_m,tip_p,t,C=tip_chord,n=nu)
    upper = numpy.array([Xu,Yu]).T
    lower = numpy.array([Xl,Yl]).T
    tip_section = numpy.concatenate((lower[::-1],upper[1:]))
    # set x position
    tip_section[:,0]-=0.25*tip_chord
    # build mesh
    vtk_model = vtk.vtkStructuredGrid()
    vtk_model.SetDimensions(2*nu-1,nv,1)
    # build points
    vtk_points = vtk.vtkPoints()
    theta = numpy.linspace(0.,numpy.pi,nv)
    y_array = -semi_span*numpy.cos(theta)
    #y_array = numpy.linspace(-semi_span,semi_span,nv)
    for j,y in enumerate(y_array):
        r = abs(y)/semi_span
        current_section = (1.-r)*root_section+r*tip_section
        for i in xrange(2*nu-1):
            vtk_points.InsertNextPoint(current_section[i,0],y,current_section[i,1])
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
    writer.SetFileName(filename)
    writer.Write()
