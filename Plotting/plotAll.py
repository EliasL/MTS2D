from makeAnimations import makeAnimations
from makePlots import makePlots
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
    collectionName = "collection.pvd"
    energyGridName = "energy_grid.csv"

    path = dataPath+subfolderName+'/'

    makePlots(path, collectionName)
    makeAnimations(path, collectionName)
    # makeEnergyField(path, energyGridName)