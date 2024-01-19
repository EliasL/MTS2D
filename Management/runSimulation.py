from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig

outPath = "/media/elias/Data/output/"

config = SimulationConfig(nx=100, ny=100, startLoad=0.15, 
                          loadIncrement=0.00001, maxLoad=2)
manager = SimulationManager(config, outPath, False)
manager.runSimulation()
manager.plot()
