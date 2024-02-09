#include <iostream>

#include "settings.h"
#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

int main(int argc, char *argv[])
{

    // GET SIMULATION CONFIGURATION
    // Check if file and output path arguments are provided
    if (argc < 2)
    {
        std::cerr << "Error! Configuration file MUST be provided!" << std::endl;
        return 1; // Return with error code
    }

    // The name of the config file is the same as the simulation name
    std::string configFile = argv[1]; // Get the file path from command line
    std::optional<std::string> dataPath;

    // The datapath is not specified in the command
    if(argc < 3){
        // Finding output path automatically
        dataPath = std::nullopt;
    } else{
        dataPath = argv[2]; // Get the output path from command line
    }
    



    // RUN THE SIMULATION

    Simulation s = Simulation(configFile, dataPath);
    s.run_simulation();
    // drawPicture();
}

