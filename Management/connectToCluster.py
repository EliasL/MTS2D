from paramiko import SSHClient, AutoAddPolicy, AuthenticationException
import subprocess
from scp import SCPClient
from icecream import ic
import os

def custom_output(*args):
    return f"ic| {' '.join(str(arg) for arg in args)}"


def uploadProject():
    try:
        # Define the rsync command as a list of arguments
        rsync_command = [
            "rsync",
            "-avz",
            "--progress",
            "--exclude", ".git",
            "--exclude", "build",
            "--exclude", "build-release",
            "/home/elias/Work/PhD/Code/1D-version1/",
            "elundheim@galois.pmmh-cluster.espci.fr:/home/elundheim/simulation/"
        ]
        
        # Run the rsync command
        subprocess.run(rsync_command, check=True)
        print("Project folder successfully uploaded.")
        
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while uploading the project: {e}")
 
uploadProject()

def connectToCluster():
    ic.configureOutput(outputFunction=custom_output)

    username = "elundheim"
    cluster_address = "galois.pmmh-cluster.espci.fr"
    key_filename = "~/Work/ssh/eliasPmmhClusterKey" 

    # Step 1: Establish an SSH connection to the cluster using Paramiko.
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(AutoAddPolicy())

    try:
        # Connect using the private key instead of a password
        ssh.connect(cluster_address, username=username, key_filename=key_filename)
        ic("SSH connection established to the cluster.")
    except AuthenticationException:
        ic("Authentication failed. Please check your SSH key.")
        exit(1)
    except Exception as e:
        ic(f"Error connecting to the cluster: {e}")
        exit(1)

    return ssh