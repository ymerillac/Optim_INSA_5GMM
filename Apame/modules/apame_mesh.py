import numpy
import vtk

class read_vtk_mesh:
    
    def __init__(self,vtk_mesh,wake_length,te_angle_detection=140.,neighbour_limit_angle=60.):
        """
        Constructor
        """
        # angle criterion to detect trailing edges
        self.__te_criterion = numpy.cos((numpy.pi/180.)*te_angle_detection)
        # angle criterion to detect adjacent panels
        self.__adjacent_panels_criterion = numpy.cos((numpy.pi/180.)*neighbour_limit_angle)
        # vtk mesh
        self.__vtk_mesh = vtk_mesh
        self.__vtk_model = None
        self.__nb_vtk_points = None
        self.__nb_vtk_cells = None
        self.__read_vtk_mesh()
        # panels
        self.__panels = None
        # wake and trailing edge
        self.__wake_length = wake_length
        self.__te_points = None
        self.__wake_points = None
        self.__wake_panels = None
        
    def __read_vtk_mesh(self):
        """
        Import vtk mesh
        """
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(self.__vtk_mesh)
        reader.Update()
        self.__vtk_model = reader.GetOutput()
        self.__nb_vtk_cells = self.__vtk_model.GetNumberOfCells()
        self.__nb_vtk_points = self.__vtk_model.GetNumberOfPoints()
        
    def __write_vtk_mesh(self):
        """
        Export to file
        """
        vtk_writer = vtk.vtkXMLPolyDataWriter()
        if vtk.VTK_MAJOR_VERSION <= 5:
            vtk_writer.SetInput(self.__vtk_model)
        else:
            vtk_writer.SetInputData(self.__vtk_model)
        vtk_writer.SetFileName(self.__vtk_mesh)
        vtk_writer.Write()
        
    def write_apame_mesh(self,fid):
        """
        Write geometry part of the apame input file
        """
        # compute panels
        self.__compute_panels()
        self.__compute_wake_panels()
        # write points
        self.__write_points(fid)
        # write panels
        self.__write_panels(fid)
        
    def __write_points(self,fid):
        """
        Write points in apame input file
        """
        # get total number of nodes
        nb_wake_nodes = len(self.__wake_points)
        nb_nodes=self.__nb_vtk_points+nb_wake_nodes
        # write in input file
        fid.write("# GEOMETRY\n# x y z                     [m]\n")
        fid.write('NODES '+str(nb_nodes)+'\n')
        for i in xrange(self.__nb_vtk_points):
            (x,y,z) = self.__vtk_model.GetPoint(i)
            fid.write(str(x)+' '+str(y)+' '+str(z)+'\n')
        # write wake nodes
        for i in xrange(nb_wake_nodes):
            x = self.__wake_points[i][0]
            y = self.__wake_points[i][1]
            z = self.__wake_points[i][2]
            fid.write(str(x)+' '+str(y)+' '+str(z)+'\n')
            
    def write_data(self,data_dict):
        """
        Store flow data in vtk mesh
        """
        for data_id,data in data_dict.items():
            nb_cases = data.shape[1]
            for i in xrange(nb_cases):
                vtk_array = vtk.vtkDoubleArray()
                vtk_array.SetName(data_id+'_'+str(i+1))
                vtk_array.SetNumberOfComponents(1)
                for value in data[:,i]:
                    vtk_array.InsertNextValue(value)
                self.__vtk_model.GetCellData().AddArray(vtk_array)
        # update vtk file
        self.__write_vtk_mesh()
                
    
    def __write_panels(self,fid):
        """
        Write panels in apame input file
        """
        # get total number of panels
        nb_panels_wo_wake = len(self.__panels)
        nb_wake_panels = len(self.__wake_panels)
        nb_panels=nb_panels_wo_wake+nb_wake_panels
        fid.write('# type node_id1 node_id2 node_id3 node_id4  elem_id1 elem_id2 elem_id3 elem_id4\n')
        fid.write('PANELS '+str(nb_panels)+'\n')
        # write panels
        for panel in self.__panels+self.__wake_panels:
            panel_type =  panel[0]
            fid.write(str(panel_type)+' ')
            pt_list = panel[1]
            for pt in pt_list:
                fid.write(str(pt+1)+' ')
            adj_cells = panel[2]
            for adj_cell in adj_cells:
                fid.write(str(adj_cell+1)+' ')
            fid.write('\n')
            
    def __compute_panels(self):
        """
        Compute panels:
        - points compsing panel
        - neighbourhood
        """
        self.__panels = []
        self.__wake_panels = []
        self.__te_points = []
        # loop over cells
        for i in xrange(self.__nb_vtk_cells):
            panel_type,pt_list,nghb_panels,te_pts = self.__compute_panel(i)
            new_panel = (panel_type,pt_list,nghb_panels)
            self.__panels.append(new_panel)
            # add wake panel if needed
            if len(te_pts)!=0:
                assert(len(te_pts)==2)
                new_wake_panel = sorted(te_pts)
                if new_wake_panel not in self.__wake_panels:
                    self.__wake_panels.append(new_wake_panel)
                for te_pt in te_pts:
                    if te_pt not in self.__te_points:
                        self.__te_points.append(te_pt)
        
    def __compute_wake_panels(self):
        """
        Compute wake panels
        """
        self.__wake_points = []
        # build wake end points
        for te_pt in self.__te_points:
            te_pt_xyz = numpy.array(self.__vtk_model.GetPoint(te_pt))
            wake_pt_xyz = numpy.copy(te_pt_xyz)
            wake_pt_xyz[0]+=self.__wake_length
            self.__wake_points.append(wake_pt_xyz)
        # complete wake panels information
        nb_wake_panels = len(self.__wake_panels)
        for i in xrange(nb_wake_panels):
            wake_pts = self.__wake_panels[i]
            p1 = wake_pts[0]
            p2 = wake_pts[1]
            i1 = self.__te_points.index(p1)
            i2 = self.__te_points.index(p2)
            p3 = i2+self.__nb_vtk_points
            p4 = i1+self.__nb_vtk_points
            wake_pts = [p1,p2,p3,p4]
            # get adjacent panels
            adj_cells = self.__get_common_cells(p1, p2)
            assert(len(adj_cells)==2)
            # get first adjacent cell normal
            ref_normal = self.__get_cell_normal(adj_cells[0])
            if not self.__check_wake_panel_orientation(wake_pts,ref_normal):
                wake_pts.reverse()
            self.__wake_panels[i] = (10,wake_pts,adj_cells)
            
    def __check_wake_panel_orientation(self,wake_pts,ref_normal):
        p1_xyz = numpy.array(self.__vtk_model.GetPoint(wake_pts[0]))
        p2_xyz = numpy.array(self.__vtk_model.GetPoint(wake_pts[1]))
        p3_xyz = self.__wake_points[wake_pts[2]-self.__nb_vtk_points]
        p4_xyz = self.__wake_points[wake_pts[3]-self.__nb_vtk_points]
        n = numpy.cross(p3_xyz-p1_xyz,p4_xyz-p2_xyz)
        return numpy.dot(n,ref_normal)>0.
        
            
    def __compute_panel(self,cell_id):
        """
        Compute Apame panel from a vtk cell
        """
        # get cell points
        cell_pt_list = self.__get_cell_pts(cell_id)
        nb_cell_pts = len(cell_pt_list)
        if nb_cell_pts==3:
            panel_type=2
        elif nb_cell_pts==4:
            panel_type=1
        else:
            raise Exception,"Error in __compute_panel: number of panel points must be equal to 3 or 4"        
        # get neighbor cells
        ngbh_cells,te_points = self.__get_adjacent_cells(cell_id,cell_pt_list)
        # complete ngbh cell list 
        while len(ngbh_cells)!=4:
            ngbh_cells.append(-1)
        return panel_type,cell_pt_list,ngbh_cells,te_points
    
    def __get_cell_pts(self,cell_id):
        """
        Get list of cell's points (Ids)
        """
        # get cell points
        cell_points = []
        cellPointIds = vtk.vtkIdList()
        self.__vtk_model.GetCellPoints(cell_id,cellPointIds)
        nb_cell_pts = cellPointIds.GetNumberOfIds()
        for j in xrange(nb_cell_pts):
            cell_points.append(cellPointIds.GetId(j))
        return cell_points
    
    def __get_cell_normal(self,cell_id):
        """
        Get normal vector of a cell
        """
        normals = self.__vtk_model.GetCellData().GetNormals()
        return normals.GetTuple(cell_id)
    
    def __get_adjacent_cells(self,cell_id,cell_pt_list):
        """
        Get all the cells adjacent to input one
        """
        ngbh_cells = []
        te_points = []
        nb_cell_pts = len(cell_pt_list)
        n1 = self.__get_cell_normal(cell_id)
        # loop over edges
        for j in xrange(nb_cell_pts):
            p1 = cell_pt_list[j]
            if(j==nb_cell_pts-1):
                p2 = cell_pt_list[0]
            else:
                p2 = cell_pt_list[j+1]
            adj_cells = self.__get_common_cells(p1, p2, exclude=cell_id)
            assert(len(adj_cells) in [0,1])
            if len(adj_cells)==1:
                adj_cell = adj_cells[0]
                n2 = self.__get_cell_normal(adj_cell)
                normals_scalprod = numpy.dot(n1,n2)
                if normals_scalprod>self.__adjacent_panels_criterion:
                    ngbh_cells.append(adj_cell)
                if normals_scalprod<self.__te_criterion:
                    te_points = [p1,p2]
        return ngbh_cells,te_points
                
    def __get_common_cells(self,p1,p2,exclude=-1):
        """
        Get all the cells containing p1 and p2
        """
        cell_list_1 = vtk.vtkIdList()
        self.__vtk_model.GetPointCells(p1,cell_list_1)
        nb_ids = cell_list_1.GetNumberOfIds()
        cell_list_1 = [cell_list_1.GetId(j) for j in xrange(nb_ids)]
        cell_list_2 = vtk.vtkIdList()
        self.__vtk_model.GetPointCells(p2,cell_list_2)
        nb_ids = cell_list_2.GetNumberOfIds()
        cell_list_2 = [cell_list_2.GetId(j) for j in xrange(nb_ids)]
        common_cells = [cell for cell in cell_list_1 if cell in cell_list_2 and cell!=exclude]
        return common_cells
                 
if __name__=='__main__':
    fid = open('naca12.inp','w')
    vtkinter = apame_mesh_writer('naca0012.vtp',10.)
    vtkinter.write_apame_mesh(fid)
    fid.close()
    
        
        
        
