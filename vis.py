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
    ren.SetBackground( 0.329412, 0.34902, 0.427451 ) #Paraview blue

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


def render(act):
    ren = vtk.vtkRenderer()
    ren.SetBackground(0.329412, 0.34902, 0.427451)  # Paraview blue
    ren.AddActor(act)

    # Create a window for the renderer of size 250x250
    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren_win.SetSize(500, 500)

    # Set an user interface interactor for the render window
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    # Start the initialization and rendering
    iren.Initialize()
    ren_win.Render()
    iren.Start()


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


def create_volume(points, num_slices, width=256, height=256):
    # https://www.vtk.org/pipermail/vtkusers/2011-February/065697.html
    # https://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/ColorCells
    poly_data = vtk.vtkPolyData()
    #poly_data.SetExtent((0, 255, 0, 255, 0, num_slices-1))
    vtkpoints = vtk.vtkPoints()  # The raw point data
    vtkverts = vtk.vtkCellArray()  # The vertices that make up the point cloud
    scalars = vtk.vtkFloatArray()  # The scalar values on each vertex
    idx = 0
    point_counter = 0
    for z in range(0, num_slices):
        for y in range(0, height):
            for x in range(0, width):
                #print(x, y, z, idx, points[idx])
                if points[point_counter] > 800:
                    vtkpoints.InsertPoint(idx, x, y, z)
                #else:
                #    vtkpoints.InsertPoint(idx, 0, 0, 0)
                    vtkverts.InsertNextCell(idx)
                    scalars.InsertNextValue(points[point_counter])
                #print(vtkverts)
                    idx += 1
                point_counter += 1
    
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

    poly_data.SetPoints(vtkpoints)
    poly_data.SetVerts(vtkverts)
    poly_data.GetPointData().SetScalars(scalars)
    mapper = vtk.vtkPolyDataMapper()

    ###MC
    mc = vtk.vtkMarchingCubes()
    mc.SetInputData(poly_data)
    mc.ComputeNormalsOn()
    mc.SetValue(0, 300)

    mapper.SetInputData(mc.GetOutput())

    mapper.ScalarVisibilityOff()

    # Take the isosurface data and create geometry
    actor = vtk.vtkLODActor()
    actor.SetNumberOfCloudPoints(1000000)
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(1, 1, 1)
    ###

    
    #mapper.SetInputData(poly_data)
    # mapper.SetLookupTable(colour_table)
    # mapper.SetScalarVisibility(1)
    # mapper.SetScalarRange(0.0, 3000)
    # mapper.SetScalarModeToUsePointData() 
    # #mapper.SetScalarRange(0, 255)
    # actor = vtk.vtkActor()
    # actor.SetMapper(mapper)
    return actor

start = 1
stop = 30
a = load_binary_values(start, stop)
act = create_volume(a, stop-start)
render(act)


def load2(img):
    # DOESN'T WORK, need the use of header file for data
    smart = vtk.vtkMetaImageReader()
    smart.SetFileName("dataset/CThead."+str(img))
    smart.SetFileDimensionality(2)
    smart.Update()

    actor = vtk.vtkImageActor()
    actor.GetMapper().SetInputConnection(smart.GetOutputPort())
    return actor

#print(load_binary_values(1))


#tiff_to_3d()
"""
# https://www.evl.uic.edu/aspale/cs526/final/3-5-2-0.htm 
voiHead = vtk.vtkExtractVOI()
voiHead.SetInputData( data )
voiHead.SetVOI( 0,255, 60,255, 0,100 )
voiHead.SetSampleRate( 1,1,1 )

iso = vtk.vtkMarchingCubes()
iso.SetInputData(voiHead.GetOutput())
iso.ComputeNormalsOn()
iso.SetValue(0, 1000)

# Create geometry
geo = vtk.vtkPolyDataMapper()
geo.SetInputData( iso.GetOutput() )
geo.ScalarVisibilityOff()

actorBone = vtk.vtkLODActor()
actorBone.SetNumberOfCloudPoints( 1000000 )
actorBone.SetMapper(geo)
actorBone.GetProperty().SetColor( 1, 1, 1 )

# Create renderer
ren = vtk.vtkRenderer()
ren.SetBackground( 0.329412, 0.34902, 0.427451 ) #Paraview blue
ren.AddActor(actorBone)

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
"""
#dim = data.GetDimensions()

# print(data)

# reader2 = vtk.vtkUnstructuredGridReader()
# reader2.SetFileName("dataset/CThead.1")
# reader2.ReadAllVectorsOn()
# reader2.ReadAllScalarsOn()
# reader2.Update()
# data2 = reader2.GetOutput()
# print(reader2.GetFileType())
# print(data2)
