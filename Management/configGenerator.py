import os

class SimulationConfig:
    def __init__(self, **kwargs):
        # Simulation Settings
        self.nx = 10  # Default = 10
        self.ny = 10  # Default = 10
        self.nrThreads = 4  # Default = 4
        self.seed = 0 # Default = 0
        self.plasticityEventThreshold = 0.2 # Default 0.2

        # Loading parameters
        self.startLoad = 0  # Default = 0.0
        self.loadIncrement = 0.01  # Default = 0.01
        self.maxLoad = 1  # Default = 1
        self.noise = 0.05 # Default = 0.05

        # Tolerances and Iterations
        self.nrCorrections = 10 # Default = 7 TODO find good default
        self.epsg = 0.0  # Default = 0.0
        self.epsf = 0.0  # Default = 0.0
        self.epsx = 0.0  # Default = 0.0
        self.maxIterations = 0  # Default = 0 (Translates to unlimited iterations)

        # Update with any provided keyword arguments
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

    def generate_name(self, withExtension=True):
        name = (
            f"S{self.nx}x{self.ny}"+
            f"L{self.startLoad},{self.loadIncrement},{self.maxLoad}"+
            f"t{self.nrThreads}n{self.noise}M{self.nrCorrections}s{self.seed}"
        )
        # Conditionally append tolerances and iterations if they are not default
        if self.epsg != 0.0:
            name += f"EpsG{self.epsg}"
        if self.epsf != 0.0:
            name += f"EpsF{self.epsf}"
        if self.epsx != 0.0:
            name += f"EpsX{self.epsx}"
        if self.maxIterations != 0:
            name += f"MaxIter{self.maxIterations}"
        if self.plasticityEventThreshold != 0.2:
            name += f"PET{self.plasticityEventThreshold}"

        if withExtension:
            # Add file extension
            name += ".conf"

        return name
    
    def get_path_and_name(self, path, withExtension=True):
        filename = self.generate_name(withExtension)
        full_path = os.path.join(path, filename)  # Corrected line
        return full_path

    def write_to_file(self, path):
        full_path=self.get_path_and_name(path)
        
        with open(full_path, 'w') as file:
            file.write("# Simulation Settings\n")
            for attr, value in self.__dict__.items():
                file.write(f"{attr} = {value} # Default = {value}\n")

        return full_path


class ConfigGenerator:
    @staticmethod
    def generateOverThreads(threads_list, **kwargs):
        return [SimulationConfig(nrThreads=threads, **kwargs) for threads in threads_list]
    @staticmethod
    def generateOverSeeds(seeds, **kwargs):
        return [SimulationConfig(seed=seed, **kwargs) for seed in seeds]
    
if __name__ == "__main__":
    import os

    conf = SimulationConfig()
    path = conf.write_to_file('build/')
    # Extract the directory part from the original path
    directory = os.path.dirname(path)
    
    # Define the new file name (keep the same directory)
    new_file_name = 'smallSimulation.conf'  # Replace 'new_filename.ext' with the new name
    
    # Construct the new path with the same directory but a new file name
    new_path = os.path.join(directory, new_file_name)
    
    # Rename the file
    os.rename(path, new_path)