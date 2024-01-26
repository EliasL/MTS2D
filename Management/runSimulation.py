from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/data2/elundheim/output/"
outPath = "/media/elias/dataStorage/output/"

config = SimulationConfig(nx=15, ny=15, startLoad=0, 
                          loadIncrement=0.05, maxLoad=2)

with SimulationManager(config, outPath, onTheCluster=False, 
                            useProfiling=False) as manager:
    manager.runSimulation()
    manager.plot()
