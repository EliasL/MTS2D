from matplotlib import pyplot as plt
import numpy as np

from settings import settings
from vtkFunctions import *

def plotEnergyOverLoad(energy, load):
    plt.plot(load, energy, label=[f"Element {i}" for i in range(len(energy[0]))])


def main(path, pvd_file):
    dataPath = path + settings["DATAFOLDERPATH"]
    vtu_files = parse_pvd_file(dataPath, pvd_file)
    S, N, E = getDataSize(dataPath, vtu_files)
    load = np.zeros((S))
    possitions = np.zeros((S, N, 3))
    stress = np.zeros((S, N, 3))
    energy = np.zeros((S, E))
   
    for i, vtu_file in enumerate(vtu_files):
        possition, stress_field, energy_field = read_vtu_data(dataPath+vtu_file)
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
    plt.savefig(path+"energy.pdf")
    #plt.show()

print("Plotting...")
# The path should be the path from work directory to the folder inside the output folder. 
main('build/output/testing/','collection.pvd')
