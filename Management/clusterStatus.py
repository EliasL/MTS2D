from connectToCluster import Servers, connectToCluster
from multiprocessing import Pool
from tabulate import tabulate

class ServerInfo:
    def __init__(self):
        self.nrTotalCores=0
        self.nrUsedCores=0
        self.nrFreeCores=0
        self.nrJobsRunning=0
        self.nrJobsWaitingInQueue=0
        self.theNodeCanAcceptMoreJobs=False

def get_server_info(ssh_client):
    # Create a ServerInfo object
    si = ServerInfo()

    # Get the total number of cores in the system
    stdin, stdout, stderr = ssh_client.exec_command("grep -c ^processor /proc/cpuinfo")
    si.nrTotalCores = int(stdout.read().decode().strip())

    # Get the total number of cores allocated to running jobs
    command_busy_cores = "squeue -t R -o '%.6C' | awk '{s+=$1} END {print s}'"
    stdin, stdout, stderr = ssh_client.exec_command(command_busy_cores)
    si.nrUsedCores = int(stdout.read().decode().strip())

    si.nrFreeCores = si.nrTotalCores-si.nrUsedCores

    # Calculate the number of jobs running and waiting in the queue
    command_jobs_running = "squeue -t R | wc -l"
    stdin, stdout, stderr = ssh_client.exec_command(command_jobs_running)
    si.nrJobsRunning = int(stdout.read().decode().strip()) - 1  # Adjust for header line
    
    command_jobs_waiting = "squeue -t PD | wc -l"
    stdin, stdout, stderr = ssh_client.exec_command(command_jobs_waiting)
    si.nrJobsWaitingInQueue = int(stdout.read().decode().strip()) - 1  # Adjust for header line
    
   # Check for exclusive job settings
    command_exclusive_jobs = "squeue -h -o '%i %t %p %C %D %R' | grep ' R ' | awk '{print $6}'"
    stdin, stdout, stderr = ssh_client.exec_command(command_exclusive_jobs)
    exclusive_job_settings = stdout.read().decode().strip().split('\n')

    # Simplified logic to set theNodeCanAcceptMoreJobs
    # Adjust this based on how you define exclusivity or constraints in your jobs
    si.theNodeCanAcceptMoreJobs = 'exclusive' not in exclusive_job_settings

    return si

def get_server_short_name(full_address):
    return full_address.split('.')[0]

# Function to add color based on the value
def colorize(value, good_value, bad_value):
    # Attempt to safely evaluate the expression to a float
    evaluated_value = eval(str(value))
    if evaluated_value is None:
        return "Invalid input"

    # Determine coloring based on comparison with good and bad values
    if good_value < bad_value:  # Lower values are better
        if evaluated_value <= good_value:
            color = "\033[92m"  # Green
        elif evaluated_value >= bad_value:
            color = "\033[91m"  # Red
        else:
            color = "\033[93m"  # Yellow
    else:  # Higher values are better
        if evaluated_value >= good_value:
            color = "\033[92m"  # Green
        elif evaluated_value <= bad_value:
            color = "\033[91m"  # Red
        else:
            color = "\033[93m"  # Yellow

    return f"{color}{value}\033[0m"

def display_server_info(servers, infos):
    # Prepare the data with all server info details
    data = []
    for server, info, nr in zip(servers, infos, range(0,len(servers))):
        server_short_name = get_server_short_name(server)
        nr_unused_cores = f"{info.nrFreeCores}/{info.nrTotalCores}"
        colored_cores = colorize(nr_unused_cores, 0.7, 0.2)
        jobs_running = info.nrJobsRunning
        jobs_waiting = colorize(info.nrJobsWaitingInQueue, 0, 2)
        #can_accept_more_jobs = "\033[92mYes\033[0m" if info.theNodeCanAcceptMoreJobs else "\033[91mNo\033[0m"
        data.append([nr, server_short_name, colored_cores, jobs_running, jobs_waiting])

    # Printing the table using tabulate
    headers = ['nr', 'Server', 'Nr Free Cores', 'Nr Jobs Running', 'Nr Jobs Waiting']
    print(tabulate(data, headers=headers, tablefmt='grid'))

def task(server):
    try:
        ssh = connectToCluster(server, False)
        info = get_server_info(ssh)
        return info
    except Exception as e:
        return f"Error: {e}"

def get_all_server_info(servers=Servers.servers):
    with Pool(processes=len(servers)) as pool:
        info = pool.map(task, servers)
        
    return servers, info

def find_server(minNrThreads):
    print("Finding available server...")
    servers, infos = get_all_server_info()

    # Filter servers with enough cores and can accept more jobs
    eligible_servers = [(server, info) for server, info in zip(servers, infos)
                        if info.nrFreeCores > minNrThreads and info.theNodeCanAcceptMoreJobs]

    if not eligible_servers:
        print("No server currently meets the core requirement and can accept more jobs.")
        # Optionally, choose the server with the fewest jobs in the queue regardless of cores
        fallback_server = min(infos, key=lambda x: x.nrJobsWaitingInQueue if x.theNodeCanAcceptMoreJobs else float('inf'))
        if fallback_server.nrJobsWaitingInQueue != float('inf'):
            server = servers[infos.index(fallback_server)]
            print("Choosing the server with the fewest jobs in the queue as a fallback.")
            return server
        else:
            print("No available server found.")
            return None

    # Sort eligible servers by fewest jobs waiting in the queue
    server, info = min(eligible_servers, key=lambda x: x[1].nrJobsWaitingInQueue)
    
    print(f"Selected {get_server_short_name(server)} with {info.nrJobsWaitingInQueue} jobs in the queue.")
    return server


if __name__ == "__main__": 
    servers, info = get_all_server_info()
    display_server_info(servers, info)