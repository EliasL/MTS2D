from simulationManager import SimulationManager
from configGenerator import ConfigGenerator, SimulationConfig
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend suitable for scripts
import matplotlib.pyplot as plt

print("Starting benchmark run...")

threads = [1,4]
threads = [1, 10, 20, 30,32,34, 40, 50, 60, 63, 64]
configs = ConfigGenerator.generateOverThreads(threads, nx=30, ny=30, startLoad=0.15, 
                          loadIncrement=0.001, maxLoad=1)

# Warmup run
manager = SimulationManager(SimulationConfig()) # Default config values
manager.runSimulation()

# Real runs
runtimes = []
for config in configs:
    manager = SimulationManager(config)
    nrRuns = 1
    times = [manager.runSimulation(build=False) for _ in range(nrRuns)]
    time = sum(times)/nrRuns
    runtimes.append(time)

plt.plot(threads, runtimes)
plt.title("Runtime over number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Time (s)")
plt.savefig("Thread effect 100x100.pdf")