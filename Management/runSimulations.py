from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig
from multiprocessing import Pool
from datetime import time

def task(config):
    try:
        manager = SimulationManager(config)
        time = manager.runSimulation(False)
        manager.plot()
    except Exception as e:
        return f"Error: {e}"
    return time


seeds = range(0,10)
configs = ConfigGenerator.generate_over_seeds(seeds, nx=10, ny=10, startLoad=0.15, 
                          loadIncrement=0.0001, maxLoad=1, threads=1)

#Build and test (Fail early)
manager = SimulationManager(SimulationConfig())
manager.runSimulation()

print(time())
with Pool(processes=len(seeds)) as pool: 
     results = pool.map(task, configs)

print(time())