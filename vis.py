import vtk

#dr = vtk.vtkDataReader()

reader = vtk.vtkTIFFReader()
#reader.SetFileName("dataset/CThead.vtk")
reader.SetFileName("dataset1/cthead-8bit001.tif")
#reader.OpenFile()
#reader.SetFileName("dataset1/cthead-8bit001.tif")
#reader.ReadAllVectorsOn()
#reader.ReadAllScalarsOn()
#reader.Update()

data = reader.GetOutput()
#depth2points = vtk.vtkDepthImageToPointCloud("dataset1/cthead-8bit001.tif")

def tiff_to_3d():  # http://vtk.1045678.n5.nabble.com/reconstruct-a-stack-of-TIFF-images-in-3D-td5719585.html
    reader = vtk.vtkTIFFReader() 
    #reader.SetFileName('dataset1/cthead-8bit001.tif') 
    #reader.Update() 

    #image = reader.GetOutput() 
    # For my test TIFF file of 640x400, I get the following 
    #image.GetExtent((0, 255, 0, 255, 0, 0))

    # So this was for a single image. For the volumetric image 
    # you now do something like the following 

    n = 20
    volume = vtk.vtkImageData() 
    volume.SetExtent((0, 255, 0, 255, 0, n-1))   # N  is the number of TIFF images 

    # Get the number of cells and points for your volume 
    ncells = volume.GetNumberOfCells() 
    npoints = volume.GetNumberOfPoints() 

    # The final piece of the puzzle is to get your TIFF data into the 
    # volume image. This is actually quite simple. What you do 
    # is create a data array to hold the data for each point 
    # in your volume image The dimension is npoints. Than you loop 
    # over each TIFF image and copy the data to the array. So something 
    # like the following: 

    array = vtk.vtkUnsignedShortArray() 
    array.SetNumberOfValues(npoints) 

    for i in range(1, n): 
        extra = "00"
        if i > 99:
            extra = ""
        elif i > 9:
            extra = "0"
        reader.SetFileName('dataset1/cthead-8bit'+extra+'%d.tif' % i) 
        reader.Update() 
        img = reader.GetOutput() 
        vals = img.GetPointData().GetArray('Tiff Scalars') 
        offset = i*n
        #for i,v in enumerate(vals): 
        #    array.SetValue(offset+i,v)
        for j in range(0, vals.GetNumberOfValues()):
            array.SetValue(offset+j, int(vals.GetTuple(j)[0]))
    # Finally we have to assign our array to the volume image. 
    volume.GetPointData().SetScalars(array)  # https://www.vtk.org/pipermail/vtkusers/2003-August/019761.html 

    print(volume.GetPointData())
    #print(array)
    # Create renderer
    ren = vtk.vtkRenderer()
    ren.SetBackground(0.329412, 0.34902, 0.427451) #Paraview blue

    #print(volume)
    #volume_mapper = vtk.vtkDataSetMapper()
    #volume_mapper.SetInputData(volume)

    g_filter = vtk.vtkStructuredGridGeometryFilter()
    g_filter.SetInputData(volume)
    volume_mapper = vtk.vtkPolyDataMapper()
    volume_mapper.SetInputConnection(g_filter.GetOutputPort())

    act = vtk.vtkActor()
    act.SetMapper(volume_mapper)
    act.GetProperty().SetPointSize(3)
    #vol = vtk.vtkVolume(volume)
    ren.AddActor(act)

    # Create a window for the renderer of size 250x250
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetSize(500, 500)

    # Set an user interface interactor for the render window
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Start the initialization and rendering
    iren.Initialize()
    renWin.Render()
    iren.Start()


class vtkCallback():
    def __init__(self):
        self.xmin = 128 # CHANGE THESE FROM HARD CODED
        self.ymin = 128
        self.zmin = 56
    def key_pressed_callback(self, obj, event):
        # https://www.cs.purdue.edu/homes/cs530/code/interactor/interactor_demo.py
        # https://www.vtk.org/Wiki/VTK/Examples/Python/Animation
        key = obj.GetKeySym()
        if key == "s":
            self.actors[0].VisibilityOff()
        elif key == 'd':
            self.actors[0].VisibilityOn()
        elif key == 'k':
            self.plane.SetNormal(1, 0, 0)
            self.xmin += 5
            self.plane.SetOrigin(self.xmin, 128, 56)
            #self.plane_points.SetPoint(0, self.xmin, 128, 56)
            #self.planes.SetPoints(self.plane_points)
            #print(self.planes.GetPoints().GetPoint(0), self.plane_points.GetPoint(0))
            #self.planes.GetPlane(0).SetOrigin(self.xmin, 0, 0)
            #self.planes.SetPlane(0, self.plane)
        elif key == 'l':
            self.plane.SetNormal(1, 0, 0)
            self.xmin -= 5
            self.plane.SetOrigin(self.xmin, 128, 56)
        elif key == 'o':
            self.plane.SetNormal(0, 1, 0)
            self.ymin += 5
            self.plane.SetOrigin(128, self.ymin, 56)
        elif key == 'p':
            self.plane.SetNormal(0, 1, 0)
            self.ymin -= 5
            self.plane.SetOrigin(128, self.ymin, 56)
        elif key == 'n':
            self.plane.SetNormal(0, 0, 1)
            self.zmin += 5
            self.plane.SetOrigin(128, 128, self.zmin)
        elif key == 'm':
            self.plane.SetNormal(0, 0, 1)
            self.zmin -= 5
            self.plane.SetOrigin(128, 128, self.zmin)
        if self.xmin < 0:
            self.xmin = 0
        if self.xmin > 256:
            self.xmin = 256
        if self.ymin < 0:
            self.ymin = 0
        if self.ymin > 256:
            self.ymin = 256
        if self.zmin < 0:
            self.zmin = 0
        if self.zmin > 113:
            self.zmin = 113

        obj.GetRenderWindow().Render()


def render(acts):
    ren = vtk.vtkRenderer()
    ren.SetBackground(0.329412, 0.34902, 0.427451)  # Paraview blue

    for act in acts:
        ren.AddActor(act)

    # Create a window for the renderer of size 500x500
    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren_win.SetSize(500, 500)

    # Set an user interface interactor for the render window
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    cb = vtkCallback()
    cb.actors = acts
    cb.renderer = ren
    iren.AddObserver("KeyPressEvent", cb.key_pressed_callback)
    # Start the initialization and rendering
    iren.Initialize()

    #acts[0].VisibilityOff()
    ren_win.Render()
    iren.Start()
    print(ren.GetNearClippingPlaneTolerance())
    


def load_binary_values(start=0, stop=0):
    # https://stackoverflow.com/questions/1035340/reading-binary-file-and-looping-over-each-byte
    # https://docs.python.org/3/library/stdtypes.html#int.from_bytes
    pixel_vals = []
    for img in range(start, stop):
        with open("dataset/CThead."+str(img), "rb") as f:
            byte = 0
            while byte != b"":
                # Do stuff with byte
                byte = f.read(2)  # Read in 2 bytes = 16 bits - ASSIGNMENT SAYS 8 BITS BUT DATA SET IS 16
                val = int.from_bytes(byte, byteorder='big', signed=False)
                pixel_vals.append(val)
        pixel_vals.pop()
    return pixel_vals


def visualise_marching_cubes_2(points, num_slices, width=256, height=256):
    # https://www.vtk.org/Wiki/VTK/Examples/Python/WriteReadVtkImageData
    # https://www.vtk.org/Wiki/VTK/Examples/Python/vtkCutter
    img_data = vtk.vtkImageData()  # Poly data to hold information about the model
    img_data.SetDimensions(width, height, num_slices)
    img_data.AllocateScalars(vtk.VTK_DOUBLE, 1)  # Need this
    point_counter = 0
    for z in range(0, num_slices):
        for y in range(0, height):
            for x in range(0, width):
                img_data.SetScalarComponentFromDouble(x, y, z, 0, points[point_counter])  # Must be double for this version of vtk (>=4.4)
                point_counter += 1

    # planes = vtk.vtkPlanes()
    # pp = [(128, 128, 56), (128, 128, 56), (128, 128, 56)]
    # plane_points = vtk.vtkPoints()
    # for i in range(0, 1):
    #     plane_points.InsertNextPoint(pp[i])
    # planes.SetPoints(plane_points)

    # np = [(1.0, 0.0, 0), (0.0, 1.0, 0), (0.0, 0.0, 1.0)]
    # normal_points = vtk.vtkDoubleArray()
    # normal_points.SetNumberOfComponents(3)
    # normal_points.SetNumberOfTuples(1)
    # for i in range(0, 1):
    #     normal_points.SetTuple3(i, np[i][0], np[i][1], np[i][2])
    # planes.SetNormals(normal_points)

    plane = vtk.vtkPlane()
    plane.SetOrigin(128, 128, 56)
    plane.SetNormal(1.0, 0.0, 0.0)

    # Skin
    mc = vtk.vtkMarchingCubes()  # Perform marching cubes
    mc.ComputeNormalsOn()
    mc.SetValue(0, 600)
    mc.SetInputData(img_data)

    clipper = vtk.vtkClipPolyData()
    clipper.SetInputConnection(mc.GetOutputPort())
    clipper.SetClipFunction(plane)

    clip_mapper = vtk.vtkPolyDataMapper()
    clip_mapper.ScalarVisibilityOff()  # Use this to make sure shows the colour given below
    clip_mapper.SetInputConnection(clipper.GetOutputPort())

    clip_actor = vtk.vtkActor()
    clip_actor.SetMapper(clip_mapper)
    clip_actor.GetProperty().SetDiffuseColor(1, 0.49, 0.25)
    clip_actor.GetProperty().SetSpecular(0.3)
    clip_actor.GetProperty().SetSpecularPower(20)
    clip_actor.GetProperty().SetOpacity(0.5)

    # Bone
    mc_bone = vtk.vtkMarchingCubes()
    mc_bone.ComputeNormalsOn()
    mc_bone.SetValue(0, 1500)
    mc_bone.SetInputData(img_data)

    clipper_bone = vtk.vtkClipPolyData()
    clipper_bone.SetInputConnection(mc_bone.GetOutputPort())
    clipper_bone.SetClipFunction(plane)

    clip_mapper_bone = vtk.vtkPolyDataMapper()
    clip_mapper_bone.ScalarVisibilityOff()
    clip_mapper_bone.SetInputConnection(clipper_bone.GetOutputPort())

    clip_actor_bone = vtk.vtkActor()
    clip_actor_bone.GetProperty().SetColor(1, 1, 0.95)
    clip_actor_bone.SetMapper(clip_mapper_bone)

    backFaces = vtk.vtkProperty()
    backFaces.SetSpecular(0.0)
    backFaces.SetDiffuse(0.0)
    backFaces.SetAmbient(1.0)
    backFaces.SetAmbientColor(0.7, 0.7, 0.65)
 
    clip_actor_bone.SetBackfaceProperty(backFaces)

    acts = [clip_actor, clip_actor_bone]
    
    # Render
    ren = vtk.vtkRenderer()
    ren.SetBackground(0.329412, 0.34902, 0.427451)  # Paraview blue

    for act in acts:
        ren.AddActor(act)

    # Create a window for the renderer of size 500x500
    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren_win.SetSize(500, 500)

    # Set an user interface interactor for the render window
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    cb = vtkCallback()
    cb.actors = acts
    cb.renderer = ren
    #cb.plane_points = plane_points
    #cb.planes = planes
    cb.plane = plane
    iren.AddObserver("KeyPressEvent", cb.key_pressed_callback)
    # Start the initialization and rendering
    iren.Initialize()

    ren_win.Render()
    iren.Start()



def visualise_marching_cubes(points, num_slices, width=256, height=256):
    poly_data = vtk.vtkPolyData()  # Poly data to hold information about the model
    hull = vtk.vtkConvexHull2D()  # Objectt o perform 2D convex hull - do this on each slice of the image
    hull_points = vtk.vtkPoints()  # The points that will be in the final hull. Combine all 2D hulls to get a 3D hull
    hull_verts = vtk.vtkCellArray()  # Ther corresponding vertices for the model
    hull_scalars = vtk.vtkFloatArray()
    point_counter = 0
    for z in range(0, num_slices):
        hull_points_slice = vtk.vtkPoints()  # THe hull points for this slice
        for y in range(0, height):
            for x in range(0, width):
                if points[point_counter] > 800 and x > 10 and x < 200 and y > 0 and y < 220:  # Some hard coded boundaries
                    hull_points_slice.InsertNextPoint(x, y, z)  # Add a point in the current slice
                point_counter += 1
        hull_points_out = vtk.vtkPoints()  # The points for the hull of the current slice
        hull.CalculateConvexHull(hull_points_slice, hull_points_out)  # Calculate the convex hull for this slice
        for i in range(0, hull_points_out.GetNumberOfPoints()):  # Add the new hull points to the overall point cloud
            point = hull_points_out.GetPoint(i)
            hull_points.InsertNextPoint(point[0], point[1], float(z))
            hull_verts.InsertNextCell(hull_points.GetNumberOfPoints())
            scalar_pos = z*height*width + y*height + x
            hull_scalars.InsertNextValue(points[0])

    mapper = vtk.vtkPolyDataMapper()

    poly_data.SetPoints(hull_points)
    poly_data.SetVerts(hull_verts)
    poly_data.GetPointData().SetScalars(hull_scalars)
    print(1)
    print(hull_points, hull_verts, poly_data.GetOutput())
    # https://www.vtk.org/Wiki/VTK/Examples/Boneyard/Cxx/PolyData/ConvexHullDelaunay3D
    #delaunay = vtk.vtkDelaunay3D()  # 3D convex hull mesh
    #delaunay.SetInputData(poly_data)
    #delaunay.Update()  # Get the 3D convex hull of the point cloud
    print(2)

    mc = vtk.vtkMarchingCubes()
    mc.ComputeNormalsOn()
    mc.SetValue(0, 20)

    volume = vtk.vtkImageData()
    volume.SetPoints(hull_points_out)
    #volume.DeepCopy(poly_data)
    print(volume)
    mc.SetInputData(volume)
    mc.ComputeNormalsOn()

    # surface = vtk.vtkSurfaceReconstructionFilter()
    # surface.SetInputData(delaunay.GetOutput())

    # cf = vtk.vtkContourFilter()
    # cf.SetInputConnection(surface.GetOutputPort())
    # cf.SetValue(0, 0.0)

    # reverse = vtk.vtkReverseSense()
    # reverse.SetInputConnection(cf.GetOutputPort())
    # reverse.ReverseCellsOn()
    # reverse.ReverseNormalsOn()

    mapper.SetInputConnection(mc.GetOutputPort())
    mapper.ScalarVisibilityOff()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(2)
    return actor


def visualise_convex_hull(points, num_slices, surface=True, width=256, height=256):
    poly_data = vtk.vtkPolyData()  # Poly data to hold information about the model
    hull = vtk.vtkConvexHull2D()  # Object to perform 2D convex hull - do this on each slice of the image
    hull_points = vtk.vtkPoints()  # The points that will be in the final hull. Combine all 2D hulls to get a 3D hull
    hull_verts = vtk.vtkCellArray()  # Ther corresponding vertices for the model
    point_counter = 0
    for z in range(0, num_slices):
        hull_points_slice = vtk.vtkPoints()  # THe hull points for this slice
        for y in range(0, height):
            for x in range(0, width):
                if points[point_counter] > 800 and x > 10 and x < 200 and y > 0 and y < 220:  # Some hard coded boundaries
                    hull_points_slice.InsertNextPoint(x, y, z)  # Add a point in the current slice
                point_counter += 1
        hull_points_out = vtk.vtkPoints()  # The points for the hull of the current slice
        hull.CalculateConvexHull(hull_points_slice, hull_points_out)  # Calculate the convex hull for this slice
        for i in range(0, hull_points_out.GetNumberOfPoints()):  # Add the new hull points to the overall point cloud
            point = hull_points_out.GetPoint(i)
            hull_points.InsertNextPoint(point[0], point[1], float(z))
            #print(hull_points.GetNumberOfPoints()-1, i)
            hull_verts.InsertNextCell(hull_points.GetNumberOfPoints()-1)

    mapper = vtk.vtkPolyDataMapper()
    #print(poly_data)
    #for i in range(0, hull_points.GetNumberOfPoints()):
    #    print(hull_points.GetPoint(i))
        #print(hull_verts.GetCell(i))
    poly_data.SetPoints(hull_points)
    #poly_data.SetVerts(hull_verts)
    #print(hull_points, hull_verts)
    #print(poly_data)

    if surface:  # If show with result with a surface instead of point cloud
        # https://www.vtk.org/Wiki/VTK/Examples/Boneyard/Cxx/PolyData/ConvexHullDelaunay3D
        delaunay = vtk.vtkDelaunay3D()  # 3D convex hull mesh
        #print(1)
        delaunay.SetInputData(poly_data)
        #print(2)
        delaunay.Update()  # Get the 3D convex hull of the point cloud
        #print(3)

        surface = vtk.vtkSurfaceReconstructionFilter()
        surface.SetInputData(delaunay.GetOutput())

        cf = vtk.vtkContourFilter()
        cf.SetInputConnection(surface.GetOutputPort())
        cf.SetValue(0, 0.0)

        reverse = vtk.vtkReverseSense()
        reverse.SetInputConnection(cf.GetOutputPort())
        reverse.ReverseCellsOn()
        reverse.ReverseNormalsOn()

        mapper.SetInputConnection(cf.GetOutputPort())
        mapper.ScalarVisibilityOff()
    else:  # Show as point cloud
        poly_data.SetVerts(hull_verts)
        mapper.SetInputData(poly_data)
        #mapper.SetLookupTable(colour_table)
        #mapper.SetScalarVisibility(1)
        #mapper.SetScalarRange(0.0, 3000)
        #mapper.SetScalarModeToUsePointData()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(2)
    return actor


def visualise_point_cloud(points, num_slices, width=256, height=256):
    # https://www.vtk.org/pipermail/vtkusers/2011-February/065697.html
    # https://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/ColorCells
    poly_data = vtk.vtkPolyData()
    vtkpoints = vtk.vtkPoints()  # The raw point data
    vtkverts = vtk.vtkCellArray()  # The vertices that make up the point cloud
    scalars = vtk.vtkFloatArray()  # The scalar values on each vertex

    idx = 0  # The index of the current point to be added to the could
    point_counter = 0  # The overall points so far - used to get the actual image value for the current coordinate
    for z in range(0, num_slices*1, 1):
        for y in range(0, height):
            for x in range(0, width):
                if points[point_counter] > 800 and x > 10 and x < 200 and y > 0 and y < 220:  # Hard coded boundaries - IMPROVE
                    vtkpoints.InsertPoint(idx, x, y, z)  # Add the point
                    vtkverts.InsertNextCell(idx)  # Associate a vertex with a point
                    scalars.InsertNextValue(points[point_counter])  # Associate with the image value
                    idx += 1
                point_counter += 1

    colour_table = vtk.vtkLookupTable()  # Create a lookup table to map point value to colour
    noc = max(points)
    colour_table.SetTableRange(0.0, noc)
    colour_table.SetNumberOfColors(noc)
    colour_table.Build()
    for i in range(0, noc):  # Assign the colours
        f = float(i)
        colour_table.SetTableValue(i, f/noc, f/noc, f/noc, 1)

    poly_data.SetPoints(vtkpoints)  # Set the data
    poly_data.SetVerts(vtkverts)
    poly_data.GetPointData().SetScalars(scalars)

    mapper = vtk.vtkPolyDataMapper()  # Create the mapper
    mapper.SetInputData(poly_data)  # Set the data
    mapper.SetLookupTable(colour_table)  # Associate the lookup table
    mapper.SetScalarVisibility(1)
    mapper.SetScalarRange(0.0, 3000)
    mapper.SetScalarModeToUsePointData()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor


def old_attempts(points, num_slices, width=256, height=256):
    """
    DEPRECATED
    """
    # https://www.vtk.org/pipermail/vtkusers/2011-February/065697.html
    # https://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/ColorCells
    poly_data = vtk.vtkPolyData()
    #poly_data.SetExtent((0, 255, 0, 255, 0, num_slices-1))
    vtkpoints = vtk.vtkPoints()  # The raw point data
    vtkverts = vtk.vtkCellArray()  # The vertices that make up the point cloud
    scalars = vtk.vtkFloatArray()  # The scalar values on each vertex
    hull = vtk.vtkConvexHull2D() ###
    hull_points = vtk.vtkPoints() ###
    hull_verts = vtk.vtkCellArray()
    idx = 0
    point_counter = 0
    for z in range(0, num_slices*1, 1):
        hull_points_slice = vtk.vtkPoints() ###
        #hull_idx = 0
        for y in range(0, height):
            for x in range(0, width):
                if points[point_counter] > 800 and x > 10 and x < 200 and y > 0 and y < 220:
                    vtkpoints.InsertPoint(idx, x, y, z)
                    vtkverts.InsertNextCell(idx)
                    scalars.InsertNextValue(points[point_counter])

                    #hull_points_slice.InsertPoint(hull_idx, x, y, z) ###
                    hull_points_slice.InsertNextPoint(x, y, z)
                    #hull_verts.InsertNextCell(hull_idx)
                    #print(x, y, z)
                    idx += 1
                    #hull_idx += 1
                point_counter += 1
        hull_points_out = vtk.vtkPoints() ###
        hull.CalculateConvexHull(hull_points_slice, hull_points_out) ###
        #print(hull_points_out.GetNumberOfPoints())
        #hull_points.InsertPoints(hull_points.GetNumberOfPoints(), hull_points_out.GetNumberOfPoints(), 0, hull_points_out) ###
        for i in range(0, hull_points_out.GetNumberOfPoints()):
            #print(i)
            point = hull_points_out.GetPoint(i)
            #print(point[0], point[1], z)
            #hull_points.InsertPoint(hull_points.GetNumberOfPoints()+i, point[0], point[1], z)
            hull_points.InsertNextPoint(point[0], point[1], float(z))
            hull_verts.InsertNextCell(hull_points.GetNumberOfPoints())
        #print("\n\n")
    #hull_verts = vtk.vtkCellArray()  ###
    # hull_scalars = vtk.vtkFloatArray()
    #for i in range(0, hull_points.GetNumberOfPoints()): ###
        #print(i)
    #    hull_verts.InsertNextCell(i) ###

    print(hull_verts, hull_points)

    colour_table = vtk.vtkLookupTable()
    # colour_table.SetNumberOfTableValues(255*5)
    noc = max(points)
    colour_table.SetTableRange(0.0, noc)
    colour_table.SetNumberOfColors(noc)
    colour_table.Build()
    for i in range(0, noc):
        f = float(i)
        # print(i, f/20.0, f/20.0, f/20.0, f)
        colour_table.SetTableValue(i, f/noc, f/noc, f/noc, 1)
    # colour_table.Build()
    #abc = [-1, -1, -1]
    #for i in range(0, noc):
        #colour_table.GetColor(i, abc)
        #print(abc)

    #poly_data.SetPoints(vtkpoints)
    #poly_data.SetVerts(vtkverts)
    #poly_data.GetPointData().SetScalars(scalars)
    mapper = vtk.vtkPolyDataMapper()

    ###MC
    if method == 0:
        # Marching cubes - doesn't work
        mc = vtk.vtkMarchingCubes()
        mc.SetValue(0, 300)
        mc.SetInputData(poly_data)
        mc.ComputeNormalsOn()
        print(mc.GetOutput(), poly_data)
        mapper.SetInputData(mc.GetOutput())

        #mapper.ScalarVisibilityOff()

        # Take the isosurface data and create geometry
        actor = vtk.vtkLODActor()
        #actor.SetNumberOfCloudPoints(1000000)
        actor.SetMapper(mapper)
        #actor.GetProperty().SetColor(1, 1, 1)
        ###
    elif method == 1:
        # Doesn't work
        surface = vtk.vtkSurfaceReconstructionFilter()
        surface.SetInputData(poly_data)

        cf = vtk.vtkContourFilter()
        cf.SetInputConnection(surface.GetOutputPort())
        cf.SetValue(0, 0.0)

        #reverse = vtk.vtkReverseSense()
        #reverse.SetInputConnection(cf.GetOutputPort())
        #reverse.ReverseCellsOn()
        #reverse.ReverseNormalsOn()

        mapper.SetInputConnection(cf.GetOutputPort())
        mapper.ScalarVisibilityOff()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
    elif method == 2:
        # 3D convex hull - slow, kinda works
        delaunay = vtk.vtkDelaunay3D()
        delaunay.SetInputData(poly_data)
        delaunay.Update()

        surface = vtk.vtkSurfaceReconstructionFilter()
        surface.SetInputData(delaunay.GetOutput())

        cf = vtk.vtkContourFilter()
        cf.SetInputConnection(surface.GetOutputPort())
        cf.SetValue(0, 0.0)

        #reverse = vtk.vtkReverseSense()
        #reverse.SetInputConnection(cf.GetOutputPort())
        #reverse.ReverseCellsOn()
        #reverse.ReverseNormalsOn()

        mapper.SetInputConnection(cf.GetOutputPort())
        mapper.ScalarVisibilityOff()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
    elif method == 3:
        # Point cloud - WORKS
        mapper.SetInputData(poly_data)
        mapper.SetLookupTable(colour_table)
        mapper.SetScalarVisibility(1)
        mapper.SetScalarRange(0.0, 3000)
        mapper.SetScalarModeToUsePointData()
        #mapper.SetScalarRange(0, 255)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
    elif method == 4:
        # 2D convex hulls
        # hull = vtk.vtkConvexHull2D()
        # hull_points = vtk.vtkPoints()
        # for i in range(0, num_slices):
        #     slice_points = vtk.vtkPoints()
        #     point_list = vtk.vtkIdList()
        #     for n in range(i*width*height, (i+1)*width*height):
        #         point_list.InsertNextId(n)
        #     print(1)
        #     vtkpoints.GetPoints(point_list, slice_points)
        #     print(2)
        #     hull_points_temp = vtk.vtkPoints()
        #     print(3)
        #     hull.CalculateConvexHull(slice_points, hull_points_temp)
        #     print(4)
        #     hull_points.InsertPoints(i*width*height, width*height, 0, hull_points_temp)
        #print(hull_points)
        poly_data.SetPoints(hull_points)
        poly_data.SetVerts(hull_verts)
        print("DONE")

        delaunay = vtk.vtkDelaunay3D()
        delaunay.SetInputData(poly_data)
        delaunay.Update()

        surface = vtk.vtkSurfaceReconstructionFilter()
        surface.SetInputData(delaunay.GetOutput())

        cf = vtk.vtkContourFilter()
        cf.SetInputConnection(surface.GetOutputPort())
        cf.SetValue(0, 0.0)

        reverse = vtk.vtkReverseSense()
        reverse.SetInputConnection(cf.GetOutputPort())
        reverse.ReverseCellsOn()
        reverse.ReverseNormalsOn()

        mapper.SetInputConnection(cf.GetOutputPort())
        mapper.ScalarVisibilityOff()

        #mapper.SetInputData(poly_data)
        #mapper.SetLookupTable(colour_table)
        #mapper.SetScalarVisibility(1)
        #mapper.SetScalarRange(0.0, 3000)
        #mapper.SetScalarModeToUsePointData()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetPointSize(2)
    return actor

start = 1
stop = 113
a = load_binary_values(start, stop+1)
#act = visualise_convex_hull(a, stop-start, True)
#act = visualise_point_cloud(a, stop-start)
act = visualise_marching_cubes_2(a, stop-start)

#render(act)
