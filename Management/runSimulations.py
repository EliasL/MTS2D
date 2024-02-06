from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig
from multiprocessing import Pool

def task(config):
    try:
        manager = SimulationManager(config)
        time = manager.runSimulation()
        manager.plot()
    except Exception as e:
        return f"Error: {e}"
    return time


seeds = range(0,10)
configs = ConfigGenerator.generate_over_seeds(seeds, nx=100, ny=100, startLoad=0.15, 
                          loadIncrement=0.00001, maxLoad=1, threads=4)

with Pool(processes=1) as pool: 
    results = pool.map(task, configs)
