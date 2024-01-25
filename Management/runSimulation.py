from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/media/elias/dataStorage/output/"
outPath = "/data2/elundheim/output/"

config = SimulationConfig(nx=50, ny=50, startLoad=0.15, 
                          loadIncrement=0.01, maxLoad=1)

with SimulationManager(config, outPath, onTheCluster=True, 
                            useProfiling=False) as manager:
    manager.runSimulation()
    manager.plot()
