from simulationManager import SimulationManager
from configGenerator import ConfigGenerator
from matplotlib import pyplot as plt
outPath = "/media/elias/T7 Sheild/output/"

threads = range(1,9)
configs = ConfigGenerator.generateOverThreads(threads, nx=100, ny=100, startLoad=0.15, 
                          loadIncrement=0.01, maxLoad=1)

# Warmup run
manager = SimulationManager(configs[-4], outPath, onTheCluster=False)
manager.runSimulation()

# Real runs
runtimes = []
for config in configs:
    manager = SimulationManager(config, outPath, onTheCluster=False)
    nrRuns = 1
    times = [manager.runSimulation() for _ in range(nrRuns)]
    time = sum(times)/nrRuns
    runtimes.append(time)

plt.plot(threads, runtimes)
plt.title("Runtime over number of threads")
plt.xlabel("Number of threads")
plt.ylabel("Time (s)")
plt.savefig("Thread effect 100x100.pdf")