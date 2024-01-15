from makeAnimations import makeAnimations
from makePlots import makePlots
from makeEnergyField import makeEnergyField
from settings import settings

if __name__ == "__main__":
    outputPath = "build-release/"+settings['OUTPUTFOLDERPATH']
    subfolderName = settings['SUBFOLDERPATH']
    collectionName = "collection.pvd"
    energyGridName = "energy_grid.csv"

    path = outputPath+subfolderName

    makePlots(path, collectionName)
    makeAnimations(path, collectionName)
    # makeEnergyField(path, energyGridName)