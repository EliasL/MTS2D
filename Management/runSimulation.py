from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/data2/elundheim/output/"
outPath = "/media/elias/Data/output/"

config = SimulationConfig(nx=50, ny=50, startLoad=0.15, 
                          loadIncrement=0.01, maxLoad=1)
manager = SimulationManager(config, outPath, onTheCluster=False, 
                            useProfiling=False)
manager.runSimulation()
manager.plot()
