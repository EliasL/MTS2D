#include <iostream>

#include "settings.h"
#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

int main(int argc, char *argv[])
{
    // GET SIMULATION CONFIGURATION
    // Check if a file path argument is provided
    if (argc < 2)
    {
        std::cerr << "Error! No configuration file provided!" << std::endl;
        return 1; // Return with error code
    }

    // The name of the config file is the same as the simulation name
    std::string configFile = argv[1]; // Get the file path from command line



    // RUN THE SIMULATION

    Simulation s = Simulation(configFile);
    s.run_simulation();

    // drawPicture();
}

