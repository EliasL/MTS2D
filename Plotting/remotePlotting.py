import os
import sys
import re
from itertools import groupby
from concurrent.futures import ThreadPoolExecutor, as_completed

# Add Management to sys.path (used to import files)
sys.path.append(str(os.Path(__file__).resolve().parent.parent / 'Management'))
# Now we can import from Management
from connectToCluster import connectToCluster, Servers, get_server_short_name



def find_csv_on_server(self, server, config):
    # Connect to the server
    ssh = connectToCluster(server, False)

    # Check if /data2 exists, otherwise use /data
    stdin, stdout, stderr = ssh.exec_command("if [ -d /data2 ]; then echo '/data2'; else echo '/data'; fi")
    base_dir = stdout.read().strip().decode()

    data_path = os.path.join(base_dir, self.user)

    # Navigate to the base_dir and get the first folder name
    command = f"cd /{data_path}; ls -d */ | head -n 1"
    stdin, stdout, stderr = ssh.exec_command(command)
    folder_name = stdout.read().strip().decode().rstrip('/')
    
    # If there is no data, we can return nothing now
    if folder_name == '':
        return []

    # Warning if the folder is not 2DCS_output
    if folder_name != "2DCS_output":
        print(f"Warning: The folder in {data_path} on {server} is not called 2DCS_output. Found: {folder_name}")

    # List all folders within the output folder
    command = f"cd /{data_path}/{folder_name}; ls -d */"
    stdin, stdout, stderr = ssh.exec_command(command)
    folders = stdout.read().strip().decode().split('\n')
    folders = [folder.rstrip('/') for folder in folders]  # Clean up folder names

    # Save the list of folders in data dictionary
    return folders

# This function searches all the servers for the given config file,
# downloads the csv file associated with the config file to a temp file,
# and returns the new local path to the csv
def get_csv_file(config):
    user="elundheim"
    newPaths = []
    # Use ThreadPoolExecutor to execute find_data_on_server in parallel across all servers
    with ThreadPoolExecutor(max_workers=len(Servers.servers)) as executor:
        future_to_server = {executor.submit(find_csv_on_server, server): server for server in Servers.servers}
        for future in as_completed(future_to_server):
            server = future_to_server[future]
            try:
                path = future.result()
                newPaths.append(path)
            except Exception as exc:
                print(f'{server} generated an exception: {exc}')

    # Process each server
    for server, folders in self.data.items():
        grouped_folders = self.parse_and_group_seeds(folders)
        if len(grouped_folders)==0:
            continue
        print(f"{get_server_short_name(server)}:")
        for folder in grouped_folders:
        print(f"\t{folder}")