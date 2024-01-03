from makeAnimations import makeAnimations
from makePlots import makePlots
from makeEnergyField import makeEnergyField

if __name__ == "__main__":
    outputPath = "build/output/"
    subfolderName = "testing/" # This name should be given by args
    collectionName = "collection.pvd"
    energyGridName = "energy_grid.csv"

    path = outputPath+subfolderName

    makePlots(path, collectionName)
    makeAnimations(path, collectionName)
    makeEnergyField(path, energyGridName)