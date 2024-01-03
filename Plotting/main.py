from makeAnimations import makeAnimations
from makePlots import makePlots

if __name__ == "__main__":
    outputPath = "build/output/"
    subfolderName = "testing/" # This name should be given by args
    collectionName = "collection.pvd"

    path = outputPath+subfolderName

    makePlots(path, collectionName)
    makeAnimations(path, collectionName)