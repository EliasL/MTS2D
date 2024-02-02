from paramiko import SSHClient, AutoAddPolicy, AuthenticationException
import subprocess
from scp import SCPClient
import os

class Servers:
    # https://intrapmmh.spip.espci.fr/spip.php?article14
    servers = [
        "galois.pmmh-cluster.espci.fr",
        "pascal.pmmh-cluster.espci.fr",
        "schwartz.pmmh-cluster.espci.fr",
        "lagrange.pmmh-cluster.espci.fr",
        "condorcet.pmmh-cluster.espci.fr",
        "dalembert.pmmh-cluster.espci.fr",
        "poincare.pmmh-cluster.espci.fr",
        "fourier.pmmh-cluster.espci.fr",
        "descartes.pmmh-cluster.espci.fr",
    ]
    default = servers[0]

def uploadProject(cluster_address=Servers.default):
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
            f"elundheim@{cluster_address}:/home/elundheim/simulation/"
        ]
        
        # Run the rsync command
        subprocess.run(rsync_command, check=True)
        print("Project folder successfully uploaded.")
        
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while uploading the project: {e}")

def connectToCluster(cluster_address=Servers.default):

    username = "elundheim"
    key_filename = "/home/elias/Work/ssh/eliasPmmhClusterKey" 

    # Step 1: Establish an SSH connection to the cluster using Paramiko.
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(AutoAddPolicy())

    try:
        # Connect using the private key instead of a password
        ssh.connect(cluster_address, username=username, key_filename=key_filename)
        print(f"SSH connection established to {cluster_address}.")
    except AuthenticationException:
        print(f"Authentication with {cluster_address} failed. Please check your SSH key.")
        exit(1)
    except Exception as e:
        print(f"Error connecting to {cluster_address}: {e}")
        exit(1)

    return ssh