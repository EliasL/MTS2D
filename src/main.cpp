#include <iostream>

#include "settings.h"
#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

int main(int argc, char *argv[])
{
    // GET SIMULATION CONFIGURATION
    // Check if file and output path arguments are provided
    if (argc < 3)
    {
        std::cerr << "Error! Configuration file and output path MUST be provided!" << std::endl;
        return 1; // Return with error code
    }

    // The name of the config file is the same as the simulation name
    std::string configFile = argv[1]; // Get the file path from command line
    std::string dataPath = argv[2]; // Get the output path from command line



    // RUN THE SIMULATION

    Simulation s = Simulation(configFile, dataPath);
    s.run_simulation();

    // drawPicture();
}

