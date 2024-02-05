from connectToCluster import connectToCluster, uploadProject, Servers
from fabric import Connection
import time
import sys
# Use select for non-blocking I/O
import select

server = Servers.servers[4]
uploadProject(server)

def run_remote_script(server_hostname, server_user, server_key_path, script_path):
    # Establish the SSH connection
    connect_kwargs = {"key_filename": server_key_path}  # Path to your SSH private key
    with Connection(host=server_hostname, user=server_user, connect_kwargs=connect_kwargs) as c:
        # Execute the remote command (your Python script)
        result = c.run(f'python3 -u {script_path}', hide=False, warn=True)

        # `hide=False` means output and errors are printed in real time
        # `warn=True` means execution won't stop on errors (similar to try/except)

        # Check the result
        if result.ok:
            print("Script executed successfully.")
        else:
            print(f"Script execution failed: {result.stderr}")

# Example usage
server_user = 'elundheim'
server_key_path = "/home/elias/.ssh/id_rsa"
script_path = '/home/elundheim/simulation/Management/benchmarking.py'

run_remote_script(server, server_user, server_key_path, script_path)