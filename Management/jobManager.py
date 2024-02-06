from runOnCluster import queue_remote_job, find_outpath_on_server
from clusterStatus import find_server
from connectToCluster import uploadProject

if __name__ == "__main__":
    minNrThreads = 11
    server = find_server(minNrThreads)
    uploadProject(server)
    script = "benchmarking.py"
    script = "runSimulations.py"
    command=f"python3 /home/elundheim/simulation/Management/{script}"

    jobId = queue_remote_job(
        server,
        command,
        "benchmk",
        minNrThreads,
    )


