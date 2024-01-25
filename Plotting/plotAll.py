from makeAnimations import makeAnimations
from makePlots import makeSinglePlot
from makeEnergyField import makeEnergyField
from settings import settings
import sys
from pathlib import Path

if __name__ == "__main__":

    if len(sys.argv) < 3:
        raise Exception("Config file and dataPath is required!")

    # We expect the argument to be path/name.conf, and we want just the name
    subfolderName = Path(sys.argv[1]).stem
    dataPath = sys.argv[2]
    
    collection = f"{settings['COLLECTIONNAME']}.pvd"
    macroData = f"{settings['MACRODATANAME']}.csv"

    path = dataPath+subfolderName+'/'

    makeSinglePlot(path+macroData)
    makeAnimations(path, collection)

    # energyGridName = "energy_grid.csv"
    # makeEnergyField(path, energyGridName)
    