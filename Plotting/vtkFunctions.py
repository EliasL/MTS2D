import xml.etree.ElementTree as ET
import vtk
from vtk.util.numpy_support import vtk_to_numpy

def parse_pvd_file(path, pvd_file):
    tree = ET.parse(path+pvd_file)
    root = tree.getroot()
    vtu_files = []

    for dataset in root.iter('DataSet'):
        vtu_files.append(dataset.attrib['file'])

    return vtu_files

def read_vtu_data(vtu_file_path):
    # Create a reader for the VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file_path)
    reader.Update()

    # Get the 'vtkUnstructuredGrid' object from the reader
    mesh = reader.GetOutput()

    # Extract Nodes
    nodes = vtk_to_numpy(mesh.GetPoints().GetData())

    # Extract Stress Field
    stress_field = vtk_to_numpy(mesh.GetPointData().GetArray("stress_field"))

    # Extract Energy Field
    energy_field = vtk_to_numpy(mesh.GetCellData().GetArray("energy_field"))



    return nodes, stress_field, energy_field


def getDataFromName(name):
    # Split the filename by underscores
    parts = name.split('_')

    # Initialize an empty dictionary
    result = {}

    # We skipp the first and last part. The first part is the 'name', the last
    # part is the type, ie .N.vtu
    result['name'] = parts[0]
    for part in parts[1:-1]:
        key, value = part.split('=')
        # Add the key-value pair to the dictionary
        result[key] = value

    return result

def getDataSize(path, vtu_files):
    """
    Calculate the size of the data from VTU files.

    Args:
        path (str): The path to the directory containing the VTU files.
        vtu_files (list of str): A list of VTU filenames.

    Returns:

        nrSteps, 
        nrNodes,
        nrElements
    """
    # Number of files
    nrSteps = len(vtu_files)

    # To find the number of nodes and elements, we need to open a file
   
    nodes, stress_field, energy_field = read_vtu_data(path+vtu_files[0]) 

    # Number of nodes
    nrNodes = len(nodes)

    # Number of elements
    nrElements = len(energy_field)

    return nrSteps, nrNodes, nrElements


def makeImages(framePath, dataPath, vtu_files, show_nodes=True, show_vectors=True, show_text=True):

    nrSteps, nrNodes, nrElements = getDataSize(dataPath, vtu_files)

    # Set up the rendering pipeline (assuming you have your data set up)
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Configure render window size
    renderWindow.SetSize(800, 600)  # Width, Height

    # Initialize the interactor and start the rendering loop for the animation
    renderWindow.Render()
    renderWindowInteractor.Initialize()

    # Create a filter to convert the render window to an image
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)

    # Set up the writer to output PNG images
    writer = vtk.vtkPNGWriter()

    # Set up a points object to store the nodes
    points = vtk.vtkPoints()

    # Set up a polydata object to store the geometry and attributes
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    # Create and configure the mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Add the actor to the renderer
    renderer.AddActor(actor)

    if show_text:
        # Create a text actor
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput("Frame: 0\nAverage Energy: 0")
        text_actor.GetTextProperty().SetFontSize(24)
        text_actor.GetTextProperty().SetColor(1.0, 1.0, 1.0)  # White color
        # Assuming 'renderWindow' is your vtkRenderWindow object
        width, height = renderWindow.GetSize()

        # Set text position to top left corner
        # Adjust the Y-coordinate by subtracting a value to account for the text height
        text_actor.SetPosition(10, height - 100) 

        # Add the text actor to the renderer
        renderer.AddActor(text_actor)

    if show_nodes:
        # Create a source for the glyphs (e.g., spheres)
        glyph_source = vtk.vtkSphereSource()
        glyph_source.SetRadius(0.1)  # Adjust the radius as needed

        # Create the glyphs
        glyph = vtk.vtkGlyph3D()
        glyph.SetSourceConnection(glyph_source.GetOutputPort())
        glyph.SetInputData(polydata)
        glyph.ScalingOff()  # Turn off scaling if your points don't have varying sizes
        glyph.Update()

        # Use the output of the glyph filter for the mapper
        glyph_mapper = vtk.vtkPolyDataMapper()
        glyph_mapper.SetInputConnection(glyph.GetOutputPort())

        glyph_actor = vtk.vtkActor()
        glyph_actor.SetMapper(glyph_mapper)

        # Add the glyph actor to the renderer instead of the simple actor
        renderer.AddActor(glyph_actor)


    if show_vectors:
        # Create a source for the arrow glyphs
        arrow_source = vtk.vtkArrowSource()

        # Create the glyphs for arrows
        arrow_glyph = vtk.vtkGlyph3D()
        arrow_glyph.SetSourceConnection(arrow_source.GetOutputPort())
        arrow_glyph.SetInputData(polydata)
        arrow_glyph.SetVectorModeToUseVector()
        arrow_glyph.SetScaleModeToScaleByVector()
        arrow_glyph.SetScaleFactor(0.1)  # Adjust scale factor as needed
        arrow_glyph.Update()

        # Use the output of the arrow glyph filter for the mapper
        arrow_glyph_mapper = vtk.vtkPolyDataMapper()
        arrow_glyph_mapper.SetInputConnection(arrow_glyph.GetOutputPort())

        arrow_glyph_actor = vtk.vtkActor()
        arrow_glyph_actor.SetMapper(arrow_glyph_mapper)

        # Add the arrow glyph actor to the renderer
        renderer.AddActor(arrow_glyph_actor)

    # Configure the camera (example configuration)
    camera = renderer.GetActiveCamera()
    camera.SetPosition(0, 0, 10)  # Adjust as needed
    camera.SetFocalPoint(0, 0, 0) # Adjust as needed

    # Configure lighting
    renderer.AutomaticLightCreationOn()

    # Configure background color
    renderer.SetBackground(0.1, 0.2, 0.3)  # Non-black color

    # Loop through the frames of the animation
    for frame, vtu_file in enumerate(vtu_files):
        # Read the data for the current frame
        nodes, stress_field, energy_field = read_vtu_data(dataPath + vtu_file)

        if show_text:
            # Calculate the average energy
            average_energy = sum(energy_field) / len(energy_field)

            d = getDataFromName(vtu_file)
            # Update the text for the current frame
            lines = [
                f"State: {d['name']}",
                f"Frame: {frame}",
                f"Load: {d['load']}",
                f"Average Energy: {average_energy:.2f}"
            ]

            # Join the lines with a newline character and set the text
            text_actor.SetInput("\n".join(lines))

        if show_nodes:        
            # Update the points object with the new node data
            points.Reset()
            points.SetNumberOfPoints(nrNodes)
            for i, node in enumerate(nodes):
                points.SetPoint(i, node)
            

            # --- Camera
            # Calculate bounds of the nodes
            points_bounds = points.GetBounds()
            
            # Optionally, calculate the center of the bounds
            center_x = (points_bounds[0] + points_bounds[1]) / 2
            center_y = (points_bounds[2] + points_bounds[3]) / 2
            center_z = (points_bounds[4] + points_bounds[5]) / 2

            # Set camera to focus on these bounds
            camera.SetFocalPoint(center_x, center_y, center_z)
            # Adjust the position of the camera as needed
            # For example, set it to look towards the center of the nodes
            # You might need to experiment with the distance to ensure all nodes are visible
            camera.SetPosition(center_x, center_y, center_z + 10)

            # Update other camera properties as needed
            # For example, setting the view up vector
            camera.SetViewUp(0, 1, 0)

        if show_vectors:
            # Update the stress field for the arrows
            stress_field_array = vtk.vtkDoubleArray()
            stress_field_array.SetNumberOfComponents(3)
            stress_field_array.SetName("stress_field")

            for stress_vector in stress_field:
                stress_field_array.InsertNextTuple(stress_vector)

            polydata.GetPointData().SetVectors(stress_field_array)


        # Make sure the mapper and render window are updated with the new data
        polydata.Modified()
        mapper.Update()

        if not show_nodes:
            # Ensure the camera is updated if necessary        
            renderer.ResetCamera()
            renderWindow.Render()

        # Update the window to image filter and write the current frame to a file
        windowToImageFilter.Modified()
        writer.SetFileName(f"{framePath}frame_{frame:04d}.png")
        writer.SetInputConnection(windowToImageFilter.GetOutputPort())
        writer.Write()

    