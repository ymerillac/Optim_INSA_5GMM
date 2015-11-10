import os
import numpy
from apame_mesh import read_vtk_mesh

def write_inp_file(out_filename,
                   mesh_interface,
                   alpha=[0.],
                   beta=[0.],
                   v=28.,
                   rho=1.225,
                   P=101325.,
                   Mach=0.,
                   origin=[0.,0.,0.],
                   wingspan=1.,
                   ref_chord=1.,
                   Sref=1.,
                   method=0,
                   farfield_dist=5.,
                   velorder=1):
    """
    Write Apame input file
    """
    fid = open(out_filename,'w')
    # write file header
    fid.write('APAME input file\nVERSION 3.1\n')
    # write airspeed, density and pressure
    fid.write('AIRSPEED '+str(v)+'\n')
    fid.write('DENSITY '+str(rho)+'\n')
    fid.write('PRESSURE '+str(P)+'\n')
    # write Mach number
    fid.write('MACH '+str(Mach)+'\n')
    # write number of cases
    nb_cases = len(alpha)
    assert(len(beta)==nb_cases)
    fid.write('CASE_NUM '+str(nb_cases)+'\n')
    # write alpha angles for each case
    for i in xrange(nb_cases):
        fid.write(str(alpha[i])+' ')
    fid.write('\n')
    # write beta angles for each case
    for i in xrange(nb_cases):
        fid.write(str(beta[i])+' ')
    fid.write('\n')
    # write geometry properties
    fid.write('WINGSPAN '+str(wingspan)+'\n')
    fid.write('MAC '+str(ref_chord)+'\n')
    fid.write('SURFACE '+str(Sref)+'\n')
    fid.write('ORIGIN *\n')
    assert(len(origin)==3)
    for i in xrange(3):
        fid.write(str(origin[i])+' ')
    fid.write('\n')
    # write solver parameters
    fid.write('METHOD '+str(method)+'\n')
    fid.write('ERROR 1e-007\n')
    fid.write('COLLDIST 1e-007\n')
    fid.write('FARFIELD '+str(farfield_dist)+'\n')
    fid.write('COLLCALC 0\n')
    fid.write('VELORDER '+str(velorder)+'\n')
    # write required outputs
    fid.write('RESULTS 1\nRES_COEF 1\nRES_FORC 1\nRES_GEOM 1\n')
    fid.write('RES_VELO 1\nRES_PRES 1\nRES_CENT 0\nRES_DOUB 1\n')
    fid.write('RES_SORC 1\nRES_VELC 0\nRES_MESH 0\nRES_STAT 0\n')
    fid.write('RES_DYNA 0\nRES_MANO 1\n')
    # write geometry
    mesh_interface.write_apame_mesh(fid)
    # end of file
    fid.close()
    
def run_case(vtk_mesh,
             wake_length,
             alpha=[0.],
             beta=[0.],
             v=28.,
             rho=1.225,
             P=101325.,
             Mach=0.,
             origin=[0.,0.,0.],
             wingspan=1.,
             ref_chord=1.,
             Sref=1.,
             method=0,
             farfield_dist=5.,
             velorder=1):
    """
    Run Apame
    """
    # read vtk mesh
    vtk_reader = read_vtk_mesh(vtk_mesh,wake_length)
    # case name
    case_name = vtk_mesh.split('.')[0]
    # write input file
    inp_filename = case_name+'.inp'
    write_inp_file(inp_filename,
                   vtk_reader,
                   alpha=alpha,
                   beta=beta,
                   v=v,
                   rho=rho,
                   P=P,
                   Mach=Mach,
                   origin=origin,
                   wingspan=wingspan,
                   ref_chord=ref_chord,
                   Sref=Sref,
                   method=method,
                   farfield_dist=farfield_dist,
                   velorder=velorder)
    # run apame
    print "*** Running Apame ..."
    os.system('apame '+case_name)
    # check run
    result_file = case_name+'.res'
    if os.path.isfile(result_file):
        print "*** Apame ran successfully !"
    else:
        print "*** Apame ran with Errors ..."
    polar,flow_data = read_out_file(result_file,None)
    numpy.savetxt(case_name+'_polar.dat',polar)
    print "*** Store output data in "+vtk_mesh
    vtk_reader.write_data(flow_data)
    print "Done !"
        
def get_line(fid,keyword):
    line = fid.readline()
    while keyword not in line:
        line = fid.readline()
    return line

def read_panel_data(fid,keyword,nb_panel,nb_cases):
    # init result array
    out_data = numpy.zeros((nb_panel,nb_cases))
    # go to ine containing keyword
    get_line(fid,keyword)
    # get data
    for i in xrange(nb_panel):
        line = fid.readline()
        values = [eval(word) for word in line.split()]
        assert(len(values)==nb_cases)
        out_data[i,:] = numpy.array(values)
    return out_data
        
def read_out_file(out_file,data_id_list):
    """
    Read apame output file and extract
    """
    # open result file
    fid = open(out_file)
    # read number of panels
    get_line(fid,'PANEL_NUM_NO_WAKE')
    line = fid.readline()
    nb_panels = eval(line)
    # read number of cases
    get_line(fid,'CASE_NUM')
    line = fid.readline()
    case_num = eval(line)
    print "Number of cases: ",case_num
    # get alpha list
    line = fid.readline()
    alpha_list = numpy.array([eval(word) for word in line.split()])
    # get beta list
    line = fid.readline()
    beta_list = numpy.array([eval(word) for word in line.split()])
    assert(len(alpha_list)==case_num and len(beta_list)==case_num)
    # convert to degrees
    rad2deg = 180./numpy.pi
    alpha_list*=rad2deg
    beta_list*=rad2deg
    # export polar
    get_line(fid,'CDRAG')
    line = fid.readline()
    cd = numpy.array([eval(word) for word in line.split()])
    get_line(fid,'CLIFT')
    line = fid.readline()
    cl = numpy.array([eval(word) for word in line.split()])
    polar = numpy.array([alpha_list,cl,cd]).T
    #return polar
    # export other data
    res_data = {}
    for data_id in ['VELOCITY','CP','DIPOLE','SOURCE']:
        res_data[data_id] = read_panel_data(fid,data_id,nb_panels,case_num)
    return polar,res_data
    
    
    
    
    
    
    