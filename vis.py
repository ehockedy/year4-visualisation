"""
gxrv46 SSAIV Visualisation Assignment
Instructions on running code are below function definitions along with example function calls
Instructions for interacing with render window also below
"""

import vtk


class vtkCallback():
    def __init__(self, x, y, z):
        self.xmin = 0  # CHANGE THESE FROM HARD CODED
        self.ymin = y/2
        self.zmin = z/2
        self.width = x
        self.height = y
        self.depth = z

    def key_pressed_callback(self, obj, event):
        # https://www.cs.purdue.edu/homes/cs530/code/interactor/interactor_demo.py
        # https://www.vtk.org/Wiki/VTK/Examples/Python/Animation
        key = obj.GetKeySym()
        if key == "s":
            self.actors[0].SetVisibility(not self.actors[0].GetVisibility())
        elif key == 'd':
            self.actors[1].SetVisibility(not self.actors[1].GetVisibility())
        elif key == 'k':
            self.plane.SetNormal(1, 0, 0)
            self.xmin += 5
            self.plane.SetOrigin(self.xmin, self.height/2, self.depth/2)
        elif key == 'l':
            self.plane.SetNormal(1, 0, 0)
            self.xmin -= 5
            self.plane.SetOrigin(self.xmin, self.height/2, self.depth/2)
        elif key == 'o':
            self.plane.SetNormal(0, 1, 0)
            self.ymin += 5
            self.plane.SetOrigin(self.width/2, self.ymin, self.depth/2)
        elif key == 'p':
            self.plane.SetNormal(0, 1, 0)
            self.ymin -= 5
            self.plane.SetOrigin(self.width/2, self.ymin, self.depth/2)
        elif key == 'n':
            self.plane.SetNormal(0, 0, 1)
            self.zmin += 5
            self.plane.SetOrigin(self.width/2, self.height/2, self.zmin)
        elif key == 'm':
            self.plane.SetNormal(0, 0, 1)
            self.zmin -= 5
            self.plane.SetOrigin(self.width/2, self.height/2, self.zmin)

        no_change = False
        if self.xmin < 0:
            self.xmin = 255
        elif self.xmin > 255:
            self.xmin = 0
        elif self.ymin < 0:
            self.ymin = 255
        elif self.ymin > 255:
            self.ymin = 0
        elif self.zmin < 0:
            self.zmin = 112
        elif self.zmin > 112:
            self.zmin = 0
        else:
            no_change = True
        if not no_change:
            self.clipper_s.SetInsideOut(not self.clipper_s.GetInsideOut())
            if self.clipper_b is not None:
                self.clipper_b.SetInsideOut(not self.clipper_b.GetInsideOut())
        # self.xmin = self.xmin % 256
        # self.ymin = self.ymin % 256
        # self.zmin = self.zmin % 113

        obj.GetRenderWindow().Render()


def load_binary_values(start=0, stop=0, path="", num_bytes=2):
    """
    Reads the binary data from the specified range of files within the given directory\n
    Arguments:\n
    start - the number of the first file to open. It will be appended to the end of the
    path parameter, and increases in value until it reaches the value in the stop parameter\n
    stop - the number of the final file to read. The function will load all files in between\n
    path - the path to the bnary data, and the constant start of each filename to be read,
    including any . for file extensions\n
    num_bytes - the number of bytes to read at a time - 1 if data is 8-bit values, 2 if data is 16 bit values   
    """
    # https://stackoverflow.com/questions/1035340/reading-binary-file-and-looping-over-each-byte
    # https://docs.python.org/3/library/stdtypes.html#int.from_bytes
    pixel_vals = []
    for img in range(start, stop):
        with open(path+str(img), "rb") as f:
            byte = 0
            while byte != b"":
                # Do stuff with byte
                byte = f.read(num_bytes)  # Read in 2 bytes = 16 bits - ASSIGNMENT SAYS 8 BITS BUT DATA SET IS 16
                val = int.from_bytes(byte, byteorder='big', signed=False)
                pixel_vals.append(val)
        pixel_vals.pop()
    return pixel_vals


def visualise_marching_cubes(points, num_slices, width=256, height=256, iso_value_1=0.1, iso_value_2=0.33, do_smoothing=False, biggest_region=True):
    """
    Main function for the assignment. Takes in the data and renders the model.\n
    Arguments:\n
    points - the list of points obtained from the binary files\n
    num_slices - number of slices i.e. binary files being processed\n
    width - the width of each input image\n
    height - height of each input image\n
    iso_value_1 - value of the first isovalue\n
    iso_value_2 - value of the second, optional isovalue. If not wanted, make 0\n
    do_smoothing - if True, smoothing will be applied to the model\n
    biggest_region - if True, will render only the biggest region of the model\n
    """
    # https://www.vtk.org/Wiki/VTK/Examples/Python/WriteReadVtkImageData
    # https://www.vtk.org/Wiki/VTK/Examples/Python/vtkCutter
    img_data = vtk.vtkImageData()  # Poly data to hold information about the model
    img_data.SetDimensions(width, height, num_slices)
    img_data.AllocateScalars(vtk.VTK_DOUBLE, 1)  # Need this
    point_counter = 0  # The current position in the 1D list of points
    max_scalar_value = 0  # Holds the largest voxel value encountered so far

    # Iterate through the data and populate a structured data array
    for z in range(0, num_slices):
        for y in range(0, height):
            for x in range(0, width):
                img_data.SetScalarComponentFromDouble(x, y, z, 0, points[point_counter])  # Must be double for this version of vtk (>=4.4)
                point_counter += 1
                if points[point_counter] > max_scalar_value:
                    max_scalar_value = points[point_counter]

    print(max_scalar_value, max_scalar_value*iso_value_1, max_scalar_value*iso_value_2)

    # Creat the plane that will be used to clip the model
    plane = vtk.vtkPlane()
    plane.SetOrigin(0, height/2, num_slices/2)
    plane.SetNormal(1.0, 0.0, 0.0)

    # Perform marching cubes
    isosurface = vtk.vtkMarchingCubes()
    isosurface.ComputeNormalsOn()
    isosurface.SetValue(0, max_scalar_value*iso_value_1)  # The first isovalue
    isosurface.SetInputData(img_data)
    to_clip = isosurface  # Need this extra variable so that we know what to put into clipper

    if biggest_region:  # Of all the disconnected parts of the model, find the biggest
        connect_filter = vtk.vtkConnectivityFilter()
        connect_filter.SetExtractionModeToLargestRegion()
        connect_filter.SetInputConnection(isosurface.GetOutputPort())
        connect_filter.SetScalarConnectivity(False)
        connect_filter.Update()
        to_clip = connect_filter

    if do_smoothing:  # Perform smoothing on the mesh
        smooth = vtk.vtkSmoothPolyDataFilter()
        smooth.SetNumberOfIterations(15)
        smooth.SetRelaxationFactor(0.1)
        smooth.FeatureEdgeSmoothingOn()
        smooth.BoundarySmoothingOn()
        smooth.SetInputConnection(to_clip.GetOutputPort())
        smooth.Update()
        to_clip = smooth

    # Clip the model to show the cross section
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputConnection(to_clip.GetOutputPort())
    clipper.SetClipFunction(plane)

    # Create a mapper to map from clipped PolyData structure to structure that can be rendered
    clip_mapper = vtk.vtkPolyDataMapper()
    clip_mapper.ScalarVisibilityOff()  # Use this to make sure shows the colour given below
    clip_mapper.SetInputConnection(clipper.GetOutputPort())

    # Add the proprty of filling in the area of the model in contact with the clipping plane
    back_faces = vtk.vtkProperty()
    back_faces.SetSpecular(0.0)
    back_faces.SetDiffuse(0.0)
    back_faces.SetAmbient(1.0)
    back_faces.SetAmbientColor(0.7, 0.7, 0.95)

    # Create the actor object - this is the thing that gets rendered
    clip_actor = vtk.vtkActor()
    clip_actor.SetMapper(clip_mapper)
    clip_actor.GetProperty().SetDiffuseColor(1, 0.49, 0.25)
    clip_actor.SetBackfaceProperty(back_faces)

    acts = [clip_actor]  # Creat the list of actors to be rendered

    # Second iso value. This may not necessarily be calculated
    # Implementation is same as first isovalue
    if iso_value_2 > 0:
        isosurface = vtk.vtkMarchingCubes()
        isosurface.ComputeNormalsOn()
        isosurface.SetValue(0, max_scalar_value*iso_value_2)
        isosurface.SetInputData(img_data)
        to_clip = isosurface

        if do_smoothing:
            smooth = vtk.vtkSmoothPolyDataFilter()
            smooth.SetInputConnection(isosurface.GetOutputPort())
            smooth.Update()
            to_clip = smooth

        clipper_2 = vtk.vtkClipPolyData()
        clipper_2.SetInputConnection(to_clip.GetOutputPort())
        clipper_2.SetClipFunction(plane)

        clip_mapper_2 = vtk.vtkPolyDataMapper()
        clip_mapper_2.ScalarVisibilityOff()
        clip_mapper_2.SetInputConnection(clipper_2.GetOutputPort())

        clip_actor_2 = vtk.vtkActor()
        clip_actor_2.GetProperty().SetColor(1, 1, 0.95)
        clip_actor_2.SetMapper(clip_mapper_2)
        clip_actor_2.SetBackfaceProperty(back_faces)

        acts.append(clip_actor_2)

    # Render
    ren = vtk.vtkRenderer()
    ren.SetBackground(0.329412, 0.34902, 0.427451)  # Paraview blue

    for act in acts:
        ren.AddActor(act)
        act.SetPosition(0, 0, 0)  # Centre the model

    # Position the camera so that it faces the model
    cam = vtk.vtkCamera()
    cam.SetPosition(-width/2, -300, 0)
    cam.SetViewUp(0, 0, -1)
    cam.SetFocalPoint(width/2, height/2, num_slices/2)
    ren.SetActiveCamera(cam)

    # Create a window for the renderer of size 500x500
    ren_win = vtk.vtkRenderWindow()
    ren_win.AddRenderer(ren)
    ren_win.SetSize(500, 500)

    # Set an user interface interactor for the render window
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)

    # Link the callback function. This allows for custom functions upon interaction withthe render window
    cb = vtkCallback(width, height, num_slices)
    cb.actors = acts
    cb.renderer = ren
    cb.clipper_s = clipper
    if iso_value_2 > 0:
        cb.clipper_b = clipper_2
    else:
        cb.clipper_b = None
    cb.plane = plane
    iren.AddObserver("KeyPressEvent", cb.key_pressed_callback)

    # Start the initialization and rendering
    iren.Initialize()
    ren_win.Render()
    iren.Start()


"""
Running instructions:
Binary data can be loaded using load_binary_data(...)
The returned value from that can be put into visualise_marching_cubes(...) which will render the model
Examples run on test data sets below
"""

"""
# MRbrain dataset
start = 1
stop = 109
a = load_binary_values(start, stop+1, "MRbrain/MRbrain.")#"dataset/CThead.")
visualise_marching_cubes(a, stop-start, 256, 256, 0.33, 0, True, True)
"""

"""
# CThead dataset
start = 1
stop = 113
a = load_binary_values(start, stop+1, "dataset/CThead.")
visualise_marching_cubes(a, stop-start, 256, 256, 0.25, 0.5, True, True)
"""

"""
# Bunny dataset
start = 1
stop = 360
a = load_binary_values(start, stop+1, "bunny/")  # no file name prefix because the binary file names are just the numbers
visualise_marching_cubes(a, stop - start, 512, 512, 0.01, 0.0, False)
"""

"""
Render window interaction instructions
- Mouse controls default behaviour such as rotation and zooming
- s toggles view of the first isovalue model on/off
- d toggles view of the second isovalue model on/off
- o & p control the clipping in the front/back direction
- k & l control the clipping in the left/right direction
- n & m control the clipping in the top/bottom direction
"""
