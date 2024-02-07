from concurrent.futures import ThreadPoolExecutor, as_completed
from connectToCluster import connectToCluster, Servers, get_server_short_name

"""
Search through all the servers and identify all the data in all the servers
"""
import os

class DataManager:
    def __init__(self) -> None:        
        self.data = {}
        self.user="elundheim"

    def find_data_on_server(self, server):
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

    def findData(self):
        # Use ThreadPoolExecutor to execute find_data_on_server in parallel across all servers
        with ThreadPoolExecutor(max_workers=len(Servers.servers)) as executor:
            future_to_server = {executor.submit(self.find_data_on_server, server): server for server in Servers.servers}
            for future in as_completed(future_to_server):
                server = future_to_server[future]
                try:
                    folders = future.result()
                    self.data[server] = folders
                except Exception as exc:
                    print(f'{server} generated an exception: {exc}')

        # Print out the data collected
        for server, folders in self.data.items():
            if folders:  # Skip servers with empty folder lists
                print(f"{get_server_short_name(server)}:")
                for folder in folders:
                    print(f"\t{folder}")  # Print each folder on its own line with a tab indentation

if __name__ == "__main__":
    dm = DataManager()
    dm.findData()