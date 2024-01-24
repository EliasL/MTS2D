from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/media/elias/T7 Sheild/output/"
outPath = "/data2/elundheim/output/"

config = SimulationConfig(nx=50, ny=50, startLoad=0.15, 
                          loadIncrement=0.01, maxLoad=1)
manager = SimulationManager(config, outPath, onTheCluster=True, 
                            useProfiling=False)
manager.runSimulation()
manager.plot()
