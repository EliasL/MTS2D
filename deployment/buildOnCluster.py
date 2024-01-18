import os
from paramiko import SSHClient, AutoAddPolicy, AuthenticationException
from scp import SCPClient
from icecream import ic

def custom_output(*args):
    return f"ic| {' '.join(str(arg) for arg in args)}"

# configure ic to use the custom output function
ic.configureOutput(outputFunction=custom_output)

# Specify the destination directory on the cluster where you want to transfer the items.
cluster_destination = "/home/elundheim/simulation/"

# Replace 'username' and 'cluster_address' with your cluster's username and address.
username = "elundheim"
cluster_address = "galois.pmmh-cluster.espci.fr"
password = os.getenv("MY_CLUSTER_PASSWORD")

if password is None:
    ic("Password environment variable not set. Please set MY_CLUSTER_PASSWORD.")
    exit(1)

# Step 1: Establish an SSH connection to the cluster using Paramiko.
ssh = SSHClient()
ssh.load_system_host_keys()
ssh.set_missing_host_key_policy(AutoAddPolicy())

try:
    ssh.connect(cluster_address, username=username, password=password)
    ic("SSH connection established to the cluster.")
except AuthenticationException:
    ic("Authentication failed. Please check your username and password.")
    exit(1)
except Exception as e:
    ic(f"Error connecting to the cluster: {e}")
    exit(1)


# Step 2: Determines files to be transfered
# Generally, this is rather fast, so we will simply overwrite all the files.
# the only exception to this is the libs/ folder which is rather large. We will
# check if this folder already exsists, and only transfer if it does not.

# Specify the source items (directories and files) on your local machine.
source_items = ["src/", "CMakeLists.txt"]
libs_path = "libs/"
lib_folder_exists = False
try:
    sftp = ssh.open_sftp()
    path_on_cluster = os.path.join(cluster_destination, libs_path)
    sftp.stat(path_on_cluster)
    lib_folder_exists = True
    ic("Skipping /libs folder. Delete folder on cluster if it needs to be updated.")
except FileNotFoundError:
    lib_folder_exists = False
if not lib_folder_exists:  
    source_items.append(libs_path)

# Convert relative source items to absolute paths.
source_items = [os.path.abspath(item) for item in source_items]

# Step 3: Use the SCPClient from the scp library to transfer items to the cluster.
try:
    with SCPClient(ssh.get_transport()) as scp:
        for src_item in source_items:
            ic(f"Transfering {src_item.split('/')[-1]}")
            scp.put(src_item, recursive=True, remote_path=os.path.join(cluster_destination, os.path.basename(src_item)))
    ic("Items transferred to the cluster.")
except Exception as e:
    ic(f"Error transferring items: {e}")
    exit(1)


# Step 4: Build
build_command = "cd simulation && mkdir -p build-release && cd build-release && cmake -DCMAKE_BUILD_TYPE=Release .. && make"

try:
    # Execute the command and capture the channel's input, output, and error streams.
    _, stdout, stderr = ssh.exec_command(build_command)
    
    # Read and print the standard output.
    output = stdout.read().decode('utf-8')
    ic("Standard Output:")
    ic(output)
    
    # Read and print the standard error.
    error = stderr.read().decode('utf-8')
    ic("Standard Error:")
    ic(error)
    
    # Check if there were any errors.
    if stderr.channel.recv_exit_status() != 0:
        ic("Build or simulation command encountered errors.")
        exit(1)
    
    ic("Build completed on the cluster.")
except Exception as e:
    ic(f"Error executing build or simulation commands on the cluster: {e}")
    exit(1)
finally:
    # Close the SSH connection when done.
    ssh.close()

ic("Script completed successfully.")
