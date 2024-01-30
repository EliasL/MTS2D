import subprocess



def uploadLibraryFolder():
    try:
        # Define the rsync command as a list of arguments
        rsync_command = [
            "rsync",
            "-avz",
            "--progress",
            "--exclude",
            "alglib",
            "/home/elias/Work/PhD/Code/1D-version1/libs/",
            "elundheim@galois.pmmh-cluster.espci.fr:/home/elundheim/simulation/CrystalSimulation/libs/"
        ]
        
        # Run the rsync command
        subprocess.run(rsync_command, check=True)
        
        print("Library folder successfully uploaded.")
        
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while uploading the library folder: {e}")

# Call the function to execute
uploadLibraryFolder()