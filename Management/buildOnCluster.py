import os
from scp import SCPClient
from icecream import ic
import threading
from connectToCluster import connectToCluster

def custom_output(*args):
    return f"ic| {' '.join(str(arg) for arg in args)}"

def read_output(stream, label):
    while True:
        line = stream.readline()
        if not line:
            break
        print(label + line.strip())


def buildOnCluster(cluster_destination, build_command):
    # configure ic to use the custom output function
    ic.configureOutput(outputFunction=custom_output)

    # Step 1: Connect to cluster
    ssh = connectToCluster()

    # Step 2: Determines files to be transfered
    # Generally, this is rather fast, so we will simply overwrite all the files.
    # the only exception to this is the libs/ folder which is rather large. We will
    # check if this folder already exsists, and only transfer if it does not.

    # Specify the source items (directories and files) on your local machine.
    source_items = ["src/", "Management/", "Plotting/", "CMakeLists.txt"]
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
            # Ensure the remote directory exists before transferring.
            # This command creates the directory and any necessary parent directories.
            ssh.exec_command(f'mkdir -p {cluster_destination}')
            for src_item in source_items:
                remote_path = os.path.join(cluster_destination, os.path.basename(src_item))
                ic(f"Transfering {src_item.split('/')[-1]}")
                scp.put(src_item, recursive=True, remote_path=remote_path)
        ic("Items transferred to the cluster.")
    except Exception as e:
        ic(f"Error transferring items: {e}")
        exit(1)


    # Step 4: Build
    try:
        # Execute the command and capture the channel's input, output, and error streams.
        _, stdout, stderr = ssh.exec_command(f"cd {cluster_destination} &&" + build_command)

        # Start two threads to read and print the standard output and standard error in real-time.
        output_thread = threading.Thread(target=read_output, args=(stdout, "Cluster: "))
        error_thread = threading.Thread(target=read_output, args=(stderr, "Error on cluster: "))

        # Wait for both threads to finish (i.e., the command execution to complete).
        output_thread.start()
        output_thread.join()
        error_thread.start()
        error_thread.join()
        
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

    ic("Build completed successfully.")


if __name__ == "__main__":
    
    build_command = "cd simulation && mkdir -p build-release && cd build-release && cmake -DCMAKE_BUILD_TYPE=Release .. && make"
    buildOnCluster(build_command)
