from runOnCluster import queue_remote_job, find_outpath_on_server
from clusterStatus import find_server
from connectToCluster import uploadProject

if __name__ == "__main__":
    minNrThreads = 65
    server = find_server(minNrThreads)
    uploadProject(server)
    outPath = find_outpath_on_server(server)
    command="python3 /simulation/Management/benchmarking.py"
    
    jobId = queue_remote_job(
        server,
        command,
        "benchmk",
        outPath,
        minNrThreads,
    )