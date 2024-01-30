from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/data2/elundheim/output/"
outPath = "/media/elias/dataStorage/output/"

config = SimulationConfig(nx=15, ny=15, startLoad=0, 
                          loadIncrement=0.05, maxLoad=2)

manager = SimulationManager(config, outPath, useProfiling=False)
manager.runSimulation()
manager.plot()
