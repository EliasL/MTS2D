import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from runOnCluster import queue_remote_job, find_outpath_on_server
from clusterStatus import find_server, Servers
from connectToCluster import uploadProject, connectToCluster
from configGenerator import SimulationConfig

class Job:
    """
    NB This does not find slurm jobs! It checks the processes running on the 
    cluster and finds all instances of CrystalSimulation running.
    """
    def __init__(self, ssh, processID, server, timeRunning) -> None:
        self.ssh = ssh
        self.name=""
        self.p_id=processID
        self.command=""
        self.server=server
        self.progres=0
        self.output_path=""
        self.configObj=None
        self.timeRunning=timeRunning

        self.getInfoFromProcess()


    def getInfoFromProcess(self):
        stdin, stdout, stderr = self.ssh.exec_command(f"ps -p {self.p_id} -o args=")
        command_line = stdout.read().decode('utf-8').strip()
        parts = command_line.split()

        self.command = command_line
        # Assuming the second and third parts of the command are what you're interested in
        config_path = parts[1]  # This seems to be new or unused; ensure it's handled as needed
        self.output_path = parts[2]  # Assuming the last part is the output path

        self.get_config_file(config_path)
        self.name=os.path.splitext(os.path.basename(config_path))[0]
        self.get_progress()

    def get_config_file(self, config_path):
        # Download the config file using SFTP
        sftp = self.ssh.open_sftp()
        local_config_filename = f"/tmp/{self.p_id}.conf" # Extract filename from path
        sftp.get(config_path, local_config_filename)  # Download the file
        sftp.close()

        # Now parse the downloaded config file
        self.configObj = SimulationConfig()
        self.configObj.parse(local_config_filename)
        os.remove(local_config_filename)

    def get_progress(self):
        remote_file_path = (self.output_path + 
                            self.name + '/' + 
                            self.name + '.log')
        # Open the remote file for reading
        # We read the second last line because the very last line might not be complete
        with self.ssh.open_sftp().file(remote_file_path, 'rb') as file:
            # Seek to the second last line
            file.seek(-2, 2)  # Seek to the second last byte
            while file.read(1) != b'\n':  # Move to the start of the last line
                file.seek(-2, 1)  # Move back one byte
            file.seek(-2, 1)  # Move back one byte
            while file.read(1) != b'\n':  # Move to the start of the second last line
                file.seek(-2, 1)  # Move back one byte

            # Read and return the second last line
            second_last_line = file.readline().decode('utf-8').strip()
            self.progres = second_last_line

    def __str__(self) -> str:
        return (f"Job: {self.name}\n"
                f"\tServer: {self.server}\n"
                #f"\tCommand: {self.command}\n"
                f"\tTime Running: {self.timeRunning}\n"
                f"\tProgress: {self.progres}\n"
                f"\tOutput Path: {self.output_path}\n"
                f"\tID: {self.p_id}\n")


class JobManager:
    def __init__(self) -> None:        
        self.jobs = []
        self.user="elundheim"

    # Function to be executed in each thread
    def find_jobs_on_server(self, server):
        local_jobs = []
        ssh = connectToCluster(server, False)
        command = f"ps -eo pid,etime,cmd | grep [C]rystalSimulation | grep -v '/bin/sh'"
        stdin, stdout, stderr = ssh.exec_command(command)
        stdout_lines = stdout.read().decode('utf-8').strip().split('\n')
        # Filter out empty lines
        stdout_lines = [line for line in stdout_lines if line.strip()]
        for line in stdout_lines:
            parts = line.split()
            p_id = parts[0]  # PID
            time_running = parts[1]  # Elapsed time
            local_jobs.append(Job(ssh, p_id, server, time_running))
        return local_jobs
    
    def findJobs(self):
        # Use ThreadPoolExecutor to execute find_jobs_on_server in parallel across all servers
        with ThreadPoolExecutor(max_workers=len(Servers.servers)) as executor:
            # Submit all servers to the executor
            future_to_server = {executor.submit(self.find_jobs_on_server, server): server for server in Servers.servers}
            for future in as_completed(future_to_server):
                server = future_to_server[future]
                try:
                    self.jobs.extend(future.result())
                except Exception as exc:
                    print(f'{server} generated an exception: {exc}')

        for job in self.jobs:
            print(job)



if __name__ == "__main__":
    minNrThreads = 4
    #server = find_server(minNrThreads)
    
    #uploadProject(server)
    script = "runSimulations.py"
    script = "benchmarking.py"
    command=f"python3 /home/elundheim/simulation/Management/{script}"

    #jobId = queue_remote_job(server, command, "benchmk", minNrThreads)
    j=JobManager()
    #j.find_jobs_on_server("lagrange.pmmh-cluster.espci.fr")
    j.findJobs()

