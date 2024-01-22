from paramiko import SSHClient, AutoAddPolicy, AuthenticationException
from scp import SCPClient
from icecream import ic
import os

def custom_output(*args):
    return f"ic| {' '.join(str(arg) for arg in args)}"

def connectToCluster():

    ic.configureOutput(outputFunction=custom_output)

    # Replace 'username' and 'cluster_address' with your cluster's username and address.
    username = "elundheim"
    cluster_address = "galois.pmmh-cluster.espci.fr"
    password = os.getenv("MY_CLUSTER_PASSWORD")

    if password is None:
        # If you don't know how to set an environment variable, run
        # 'export MY_CLUSTER_PASSWORD=' 
        # followed by your password in a terminal
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
    return ssh