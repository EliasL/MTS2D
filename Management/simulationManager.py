from buildOnCluster import buildOnCluster
from connectToCluster import connectToCluster
import threading
import subprocess
import time
from icecream import ic

# def custom_output(*args):
#     return f"ic| {' '.join(str(arg) for arg in args)}"

# # configure ic to use the custom output function
# ic.configureOutput(outputFunction=custom_output)

class SimulationManager:

    """
    Okay, so this is a bit hard to keep track of... but here is an overview:
    Regardless of whether or not we want to run on the local machine or the
    cluster, we want to run cmake localy, upload the build folder, and run the
    final make command on the cluster (if we should run on the cluster). The 
    following tries to keep track of all the possible directories:
        (local/cluster)/(debug/release) --> outputPath    
    """

    def __init__(self, configObj, outputPath, onTheCluster, debugBuild=False, useProfiling=False):
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
        self.debugBuild = debugBuild
        self.release_build_folder = "build-release/"
        self.profile_build_folder = "build-debug/"
        self.build_folder = self.profile_build_folder if debugBuild else self.release_build_folder
        # pre-Build path (Always local)
        self.prebuild_path = self.local_destination + self.build_folder
        # Build path
        self.build_path = self.project_path + self.build_folder
        # Program path
        self.program_path = self.build_path + "CrystalSimulation"
         
        # This is set by the runSimulation method
        self.conf_file = None

    def __enter__(self):        
        # Connection to the cluster
        self.ssh = connectToCluster() if self.onTheCluster else None
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.onTheCluster:
            self.ssh.close()

    def runSimulation(self):
        self._build()
        # This function writes to the cluster if self.ssh is not None
        self.conf_file = self.configObj.write_to_file(self.build_path, self.ssh)
        run_command = f"{self.program_path} {self.conf_file} {self.outputPath}"
        if self.useProfiling:
            run_command = "valgrind --tool=callgrind " + run_command

        # Start the timer right before running the command
        start_time = time.time()
        print("Running simulation")
        self._run_command(run_command)

        # Stop the timer right after the command completes
        end_time = time.time()

        # Calculate the duration
        duration = end_time - start_time
        return duration

    def _build(self):
        print("Building...")
        build_type = "Release"#"Debug" if self.useProfiling else "Release"
        prep_command = f"mkdir -p {self.build_folder} && cd {self.build_folder} && cmake -DCMAKE_BUILD_TYPE={build_type} .. "
        build_command = prep_command + "&& make"
        if self.onTheCluster: 
            # We run a local prebuild so that the build folder is ready
            # to be uploaded to the cluster
            self._run_command(prep_command, False)
            buildOnCluster(self.cluster_destination, self.build_folder, self.ssh)
        else:
            self._run_command(build_command)
        
        print("Build completed successfully.")
    
    def plot(self):

        plotCommand = f"python3 {self.project_path}Plotting/plotAll.py {self.conf_file} {self.outputPath}"
        self._run_command(plotCommand)

    def _run_command(self, command, runLocal=None):
        # If runLocal is set, we only care about that variable
        if runLocal is not None:
            oldState = self.onTheCluster
            self.onTheCluster = runLocal
            self._run_command(command)
            self.onTheCluster = oldState
            return
        
        if not self.onTheCluster:
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
        else:
            try:

                # Execute the command and capture the channel's input, output, and error streams.
                _, stdout, stderr = self.ssh.exec_command(command)

                # Start two threads to read and print the standard output and standard error in real-time.
                output_thread = threading.Thread(target=read_output, args=(stdout, "Cluster: "))
                error_thread = threading.Thread(target=read_output, args=(stderr, "Error on cluster: "))

                # Wait for both threads to finish (i.e., the command execution to complete).
                output_thread.start()
                output_thread.join()
                error_thread.start()
                error_thread.join()
                
                # Check if there were any errors.
                if stderr.channel.recv_exit_status() != 0:
                    ic("Build or simulation command encountered errors.")
                    exit(1)
                
                ic("Build completed on the cluster.")
            except Exception as e:
                ic(f"Error executing build or simulation commands on the cluster: {e}")
                exit(1)

def read_output(stream, label):
    while True:
        line = stream.readline()
        if not line:
            break
        print(label + line.strip())
