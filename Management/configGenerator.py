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
            name += f"_EpsG{self.epsg}"
        if self.epsf != 0.0:
            name += f"_EpsF{self.epsf}"
        if self.epsx != 0.0:
            name += f"_EpsX{self.epsx}"
        if self.maxIterations != 0:
            name += f"_MaxIter{self.maxIterations}"
        if self.plasticityEventThreshold != 0.05:
            name += f"_PET{self.plasticityEventThreshold}"

        if withExtension:
            # Add file extension
            name += ".conf"

        return name

    def write_to_file(self, path):
        filename = self.generate_name()
        with open(path+filename, 'w') as file:
            file.write("# Simulation Settings\n")
            # Iterate over each attribute and its value
            for attr, value in self.__dict__.items():
                file.write(f"{attr} = {value} # Default = {value}\n")
        return path+filename


class ConfigGenerator:
    @staticmethod
    def generateOverThreads(threads_list, **kwargs):
        return [SimulationConfig(nrThreads=threads, **kwargs) for threads in threads_list]
    @staticmethod
    def generateOverSeeds(seeds, **kwargs):
        return [SimulationConfig(seed=seed, **kwargs) for seed in seeds]
    
if __name__ == "__main__":
    conf = SimulationConfig()
    conf.write_to_file('build-debug/')