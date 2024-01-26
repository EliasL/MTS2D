import os
from scp import SCPClient
import threading
from pathlib import Path

def custom_output(*args):
    return f"ic| {' '.join(str(arg) for arg in args)}"

def read_output(stream, label):
    while True:
        line = stream.readline()
        if not line:
            break
        print(label + line.strip())


def buildOnCluster(cluster_destination, build_folder, ssh):

    # Specify the source items (directories and files) on your local machine.
    source_items = [
        "Management/", 
        "Plotting/",
        build_folder
    ]
    print("Uploading files...")
    work_dir = Path(__file__).parent.parent.absolute()
    # Construct absolute paths by appending each item to the script's directory
    source_items_absolute = [str(work_dir / item) for item in source_items]
    # Step 3: Use the SCPClient from the scp library to transfer items to the cluster.
    try:

        with SCPClient(ssh.get_transport()) as scp:
            # Ensure the remote directory exists before transferring.
            # This command creates the directory and any necessary parent directories.
            ssh.exec_command(f'mkdir -p {cluster_destination}')
            for src_item in source_items_absolute:
                remote_path = os.path.join(cluster_destination, os.path.basename(src_item))
                print(f"Transfering {src_item.split('/')[-1]}")
                scp.put(src_item, recursive=True, remote_path=remote_path)
        print("Items transferred to the cluster.")
    except Exception as e:
        print(f"Error transferring items: {e}")
        exit(1)

    print("Building files...")
    # Step 4: Build
    try:

        # Execute the command and capture the channel's input, output, and error streams.
        remote_build_path = os.path.join(cluster_destination, build_folder)
        _, stdout, stderr = ssh.exec_command(f"cd {remote_build_path} && make")

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
            print("Build or simulation command encountered errors.")
            exit(1)
        
        print("Build completed on the cluster.")
    except Exception as e:
        print(f"Error executing build or simulation commands on the cluster: {e}")
        exit(1)



if __name__ == "__main__":
    
    build_command = "cd simulation && mkdir -p build-release && cd build-release && cmake -DCMAKE_BUILD_TYPE=Release .. && make"
    buildOnCluster(build_command)
