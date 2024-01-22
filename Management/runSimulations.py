from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/media/elias/Data/output/"

seeds = range(0,11)
configs = ConfigGenerator.generateOverSeeds(seeds, nx=100, ny=100, startLoad=0.15, 
                          loadIncrement=0.00001, maxLoad=1)

for config in configs:
    manager = SimulationManager(config, outPath, onTheCluster=False)
    manager.runSimulation()
    manager.plot()
