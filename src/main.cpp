#include <iostream>

// Log settings
//#define ELPP_DISABLE_LOGS
//#define ELPP_DISABLE_DEBUG_LOGS
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

int main(){ 
    
    // We fix the random seed to get reproducable results
    srand(0);
    run_simulation();
    // drawPicture();
}