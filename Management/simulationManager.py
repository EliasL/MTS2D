from pathlib import Path
import subprocess
import time
import os
import shlex

class SimulationManager:
    
    def __init__(self, configObj, outputPath=None, debugBuild=False, useProfiling=False):
        self.configObj = configObj
        self.outputPath = findOutputPath() if outputPath is None else outputPath

        self.useProfiling = useProfiling        
        self.project_path = str(Path(__file__).resolve().parent.parent)

        # Build folder
        self.debugBuild = debugBuild
        self.release_build_folder = "build-release/"
        self.profile_build_folder = "build/"
        self.build_folder = self.profile_build_folder if debugBuild else self.release_build_folder
        # Build path
        self.build_path = os.path.join(self.project_path, self.build_folder)
        # Program path
        self.program_path = self.build_path + "CrystalSimulation"
         
        # This is set by the runSimulation method
        self.conf_file = None

        # Move to project# Change the working directory
        os.chdir(self.project_path)


    def runSimulation(self):
        self._build()
        self.conf_file = self.configObj.write_to_file(self.build_path)
        simulation_command = f"{self.program_path} {self.conf_file} {self.outputPath}"
        if self.useProfiling:
            simulation_command = "valgrind --tool=callgrind " + simulation_command

        # Start the timer right before running the command
        start_time = time.time()
        print("Running simulation")
        run_command(simulation_command)

        # Stop the timer right after the command completes
        end_time = time.time()

        # Calculate the duration
        duration = end_time - start_time
        return duration


    def _build(self):
        print("Building...")
        build_type = "Release"#"Debug" if self.useProfiling else "Release"
        build_command = f"mkdir -p {self.build_folder} && cd {self.build_folder} && cmake -DCMAKE_BUILD_TYPE={build_type} .. && make"

        run_command(build_command)
        print("Build completed successfully.")
    

    def plot(self):
        plot_script = os.path.join(self.project_path, "Plotting/plotAll.py")
        plot_command = f"python3 {plot_script} {self.conf_file} {self.outputPath}"
        run_command(plot_command)


def run_command(command):
    process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline().rstrip().decode('utf-8')
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc

def findOutputPath():
    # Define the paths to check
    paths = ["/media/elias/dataStorage/output/", "/data2/elundheim/output/"]

    # Initialize a variable to store the chosen path
    chosen_path = None

    # Iterate through the paths and check if they exist
    for path in paths:
        if os.path.exists(path):
            chosen_path = path
            break  # Stop the loop once a valid path is found

    # Check if a valid path was found or raise an error
    if chosen_path is None:
        raise FileNotFoundError("None of the provided paths exist.")
    else:
        print(f"Chosen output path: {chosen_path}")
    return chosen_path