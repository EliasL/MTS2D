import re
from itertools import groupby
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
        fullPath=f"{data_path}/{folder_name}"
        command = f"cd {fullPath}; ls -d */"
        stdin, stdout, stderr = ssh.exec_command(command)
        folders = stdout.read().strip().decode().split('\n')
        folders = [os.path.join(fullPath, folder.rstrip('/')) for folder in folders]  # Clean up folder names

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

    def printData(self):
        # Process each server
        for server, folders in self.data.items():
            grouped_folders = self.parse_and_group_seeds(folders)
            if len(grouped_folders)==0:
                continue
            print(f"{get_server_short_name(server)}:")
            for folder in grouped_folders:
                print(f"\t{folder}")
    
    def delete_all_found_data(self, dryRun=True):
        for server, folders in self.data.items():
           if folders:
            self.delete_data_on_server(server, folders, dryRun)

    def delete_data_on_server(self, server, folders, dryRun=True):
        # Connect to the server
        ssh = connectToCluster(server, False)

        print(f"Are you sure you want to delete these folders on {server}?:") 
        [print(folder) for folder in folders]
        if input(f"yes/no : ")!="yes":
            return
        
        # Deleting each folder found by find_data_on_server
        for folder in folders:

            # Command to delete the folder
            delete_command = f"rm -r {folder}" 


            # Execute the delete command
            if not dryRun:
                stdin, stdout, stderr = ssh.exec_command(delete_command)
                # Check for errors
                errors = stderr.read().decode().strip()
            else:
                errors = None

            if errors:
                print(f"Error deleting {folder} on {server}: {errors}")
            else:
                print(f"Successfully deleted {folder} on {server}")

    def parse_and_group_seeds(self, folders):

        # Regex to match the base part of the folder name and the seed
        pattern = re.compile(r'(.*)s(\d+)$')

        # Parse folder names into base names and seeds
        parsed_folders = []
        for folderPath in folders:
            folder = folderPath.split('/')[-1]
            match = pattern.match(folder)
            if match:
                base_name, seed = match.groups()
                parsed_folders.append((base_name, int(seed), folder))

        # Sort by base name and seed to ensure correct grouping
        parsed_folders.sort(key=lambda x: (x[0], x[1]))

        grouped_folders = {}
        # Group by base name
        for base_name, group in groupby(parsed_folders, key=lambda x: x[0]):
            # Extract and sort seeds within each base name group
            seeds = [item for item in group]
            # Group consecutive seeds
            grouped_seeds = []
            for k, g in groupby(enumerate(seeds), lambda i_x: i_x[0] - i_x[1][1]):
                seq = list(g)
                if len(seq) > 1:
                    start, end = seq[0][1], seq[-1][1]
                    grouped_seeds.append(f"{start[1]}-{end[1]}")
                else:
                    single = seq[0][1]
                    grouped_seeds.append(str(single[1]))
            # Reconstruct the folder names for each base name group
            for seed_group in grouped_seeds:
                original_folder = seeds[0][2].rsplit('s', 1)[0]  # Get one example folder and remove seed
                grouped_folder = f"{original_folder}s{seed_group}"
                if base_name in grouped_folders:
                    grouped_folders[base_name].append(grouped_folder)
                else:
                    grouped_folders[base_name] = [grouped_folder]

        # Flatten the grouped_folders dictionary to a list
        final_grouped_folders = []
        for base_name, folders in grouped_folders.items():
            final_grouped_folders.extend(folders)

        return final_grouped_folders
    


if __name__ == "__main__":
    dm = DataManager()
    dm.findData()
    dm.printData()
    #dm.delete_all_found_data(dryRun=False)