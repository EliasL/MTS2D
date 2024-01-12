from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from icecream import ic
import os

from settings import settings
from vtkFunctions import *

def plotEnergyOverLoad(energy, load):
    plt.plot(load, energy)


def makePlots(path, pvd_file):
    ic("Plotting...")

    dataPath = path + settings["DATAFOLDERPATH"]
    if(not os.path.exists(dataPath+pvd_file)):
        print(f"No file found at: {dataPath+pvd_file}")
        return
    
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

    energy = energy.sum(axis=1) / len(energy[0])

    plotEnergyOverLoad(energy, load)

    plt.xlabel(r'$\alpha$')
    plt.ylabel('Energy')
    plt.title(r'Average energy over stress $\alpha$')

    # Automatically adjust the y-axis label position
    ax = plt.gca()
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # Optional: Makes y-ticks integers
    #ax.relim()
    ax.autoscale_view()
    plt.savefig(path+"energy.pdf")
    #plt.show()

if __name__ == "__main__":
    # The path should be the path from work directory to the folder inside the output folder. 
    makePlots('build/output/testing/','collection.pvd')
