import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

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

def plotEnergyOverLoad(energy, load):
    plt.plot(load, energy, label=[f"Element {i}" for i in range(len(energy[0]))])

def getDataFromName(name):
    # Split the filename by underscores
    parts = name.split('_')

    # Initialize an empty dictionary
    result = {}

    # We skipp the first and last part. The first part is the 'name', the last
    # part is the type, ie .N.vtu
    for part in parts[1:-1]:
        key, value = part.split('=')
        # Add the key-value pair to the dictionary
        result[key] = value

    return result

def getDataSize(path, vtu_files):
    # Number of files
    nrSteps = len(vtu_files)

    # To find the number of nodes and elements, we need to open a file
   
    nodes, stress_field, energy_field = read_vtu_data(path+vtu_files[0]) 

    # Number of nodes
    nrNodes = len(nodes)

    # Number of elements
    nrElements = len(energy_field)

    return nrSteps, nrNodes, nrElements

def main(path, pvd_file):
    vtu_files = parse_pvd_file(path, pvd_file)
    S, N, E = getDataSize(path, vtu_files)
    load = np.zeros((S))
    possitions = np.zeros((S, N, 3))
    stress = np.zeros((S, N, 3))
    energy = np.zeros((S, E))
    
    for i, vtu_file in enumerate(vtu_files):
        possition, stress_field, energy_field = read_vtu_data(path+vtu_file)
        dictData = getDataFromName(vtu_file)
        load[i] = dictData["load"]
        possitions[i] = possition
        stress[i] = stress_field
        energy[i] = energy_field

    plotEnergyOverLoad(energy, load)

    plt.xlabel(r'$\alpha$')
    plt.ylabel('Energy')
    plt.yscale('log')
    plt.title(r'Energy over stress $\alpha$')
    plt.legend()
    plt.show()

# Replace 'your_pvd_file.pvd' with the path to your .pvd file
main('../../build/output/testing/','collection.pvd')
