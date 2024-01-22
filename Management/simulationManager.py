from buildOnCluster import buildOnCluster
import subprocess
import time
from icecream import ic

def custom_output(*args):
    return f"ic| {' '.join(str(arg) for arg in args)}"

# configure ic to use the custom output function
ic.configureOutput(outputFunction=custom_output)

class SimulationManager:

    def __init__(self, configObj, outputPath, onTheCluster, useProfiling=False):
        self.configObj = configObj
        self.outputPath = outputPath
        self.onTheCluster = onTheCluster

        self.useProfiling = useProfiling        

        # Specify the destination directory on the cluster where you want to transfer the items.
        self.cluster_destination = "/home/elundheim/simulation/"
        # Same but for personal/local computer
        self.local_destination = "/home/elias/Work/PhD/Code/1D-version1/"
        self.project_path = self.cluster_destination if onTheCluster else self.local_destination
        # Build folder
        self.release_build_folder = "build-release/"
        self.profile_build_folder = "build-debug/"
        self.build_folder = self.profile_build_folder if False else self.release_build_folder
        # Build path
        self.build_path = self.project_path + self.build_folder
        # Program path
        self.program_path = self.build_path + "CrystalSimulation"
         
        # This is set by the runSimulation method
        self.conf_file = None


    def runSimulation(self):
        self._build()
        self.conf_file = self.configObj.write_to_file(self.build_path)
        run_command = f"{self.program_path} {self.conf_file} {self.outputPath}"
        if self.useProfiling:
            run_command = "valgrind --tool=callgrind " + run_command

        # Start the timer right before running the command
        start_time = time.time()

        self._run_command(run_command)

        # Stop the timer right after the command completes
        end_time = time.time()

        # Calculate the duration
        duration = end_time - start_time
        return duration

    def _build(self):
        ic("Building...")
        build_type = "Release"#"Debug" if self.useProfiling else "Release"
        build_command = f"mkdir -p {self.build_folder} && cd {self.build_folder} && cmake -DCMAKE_BUILD_TYPE={build_type} .. && make"
        if self.onTheCluster:
            buildOnCluster(self.cluster_destination, build_command)
        else:
            self._run_command(build_command)
    
    def plot(self):
        plotCommand = f"python {self.project_path}Plotting/plotAll.py {self.conf_file} {self.outputPath}"
        self._run_command(plotCommand)

    @staticmethod
    def _run_command(command):
        # Use subprocess.run to execute the command.
        result = subprocess.run(command, shell=True, text=True)

        # Get the standard output and error.
        output = result.stdout
        error = result.stderr

        # Check if the command was executed successfully
        if result.returncode == 0:
            print("Command executed successfully!")
            print("Output:\n", output)
        else:
            print("Error in command execution.")
            raise Exception("Error:\n" + error)