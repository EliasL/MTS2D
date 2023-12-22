from makeAnimations import makeAnimations
from makePlots import makePlots

if __name__ == "__main__":
    outputPath = "build/output/"
    subfolderName = "testing/"
    collectionName = "collection.pvd"

    path = outputPath+subfolderName

    makePlots(path, collectionName)
    makeAnimations(path, collectionName)