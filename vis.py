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
    print(1)
    # Start the initialization and rendering
    print(2)
    iren.Initialize()
    print(3)
    ren_win.Render()
    print(4)
    iren.Start()
    print(5)


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
    img_data = vtk.vtkImageData()  # Poly data to hold information about the model
    img_data.SetDimensions(width, height, num_slices)
    img_data.AllocateScalars(vtk.VTK_DOUBLE, 1)
    point_counter = 0
    for z in range(0, num_slices):
        for y in range(0, height):
            for x in range(0, width):
                img_data.SetScalarComponentFromDouble(x, y, z, 0, points[point_counter])  # Must be double for this version of vtk (>=4.4)
                #print(points[point_counter], img_data.GetScalarComponentAsDouble(x, y, z, 0))
                point_counter+=1

    mapper = vtk.vtkPolyDataMapper()

    mc = vtk.vtkMarchingCubes()
    mc.ComputeNormalsOn()
    mc.SetValue(0, 800)
    mc.SetInputData(img_data)
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

render([act])
