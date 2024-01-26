from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/media/elias/dataStorage/output/"

seeds = range(0,1)
configs = ConfigGenerator.generateOverSeeds(seeds, nx=100, ny=100, startLoad=0.15, 
                          loadIncrement=0.00001, maxLoad=1)

for config in configs:
    manager = SimulationManager(config, outPath, onTheCluster=True)
    manager.runSimulation()
    manager.plot()
