from pathlib import Path
import subprocess
import time


class SimulationManager:


    def __init__(self, configObj, outputPath, debugBuild=False, useProfiling=False):
        self.configObj = configObj
        self.outputPath = outputPath

        self.useProfiling = useProfiling        
        self.project_path = str(Path(__file__).resolve().parent.parent)

        # Build folder
        self.debugBuild = debugBuild
        self.release_build_folder = "build-release/"
        self.profile_build_folder = "build-debug/"
        self.build_folder = self.profile_build_folder if debugBuild else self.release_build_folder
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
        build_command = f"mkdir -p {self.build_folder} && cd {self.build_folder} && cmake -DCMAKE_BUILD_TYPE={build_type} .. && make"

        self._run_command(build_command)
        print("Build completed successfully.")
    

    def plot(self):
        plotCommand = f"python3 {self.project_path}Plotting/plotAll.py {self.conf_file} {self.outputPath}"
        self._run_command(plotCommand)


    def _run_command(self, command):
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