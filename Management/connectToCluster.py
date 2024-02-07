from paramiko import SSHClient, AutoAddPolicy, AuthenticationException
import subprocess

class Servers:
    # https://intrapmmh.spip.espci.fr/spip.php?article14
    servers = [
        "galois.pmmh-cluster.espci.fr",        #0 
        "pascal.pmmh-cluster.espci.fr",        #1 
        "schwartz.pmmh-cluster.espci.fr",      #2
        "lagrange.pmmh-cluster.espci.fr",      #3
        "condorcet.pmmh-cluster.espci.fr",     #4
        "dalembert.pmmh-cluster.espci.fr",     #5
        "poincare.pmmh-cluster.espci.fr",      #6
        "fourier.pmmh-cluster.espci.fr",       #7
        "descartes.pmmh-cluster.espci.fr",     #8
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

def connectToCluster(cluster_address=Servers.default, verbose=True):

    username = "elundheim"
    key_filename = "/home/elias/.ssh/id_rsa" 

    # Step 1: Establish an SSH connection to the cluster using Paramiko.
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(AutoAddPolicy())

    try:
        # Connect using the private key instead of a password
        ssh.connect(cluster_address, username=username, key_filename=key_filename)
        if verbose:
            print(f"SSH connection established to {cluster_address}.")
    except AuthenticationException:
        raise AuthenticationFailedException(f"Authentication with {cluster_address} failed. Please check your SSH key.")
    except Exception as e:
        raise SSHConnectionException(f"Error connecting to {cluster_address}: {e}")


    return ssh

class AuthenticationFailedException(Exception):
    pass

class SSHConnectionException(Exception):
    pass
