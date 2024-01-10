#include <iostream>

// Log settings
//#define ELPP_DISABLE_LOGS
//#define ELPP_DISABLE_DEBUG_LOGS
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

int main(){ 
    run_simulation();
    drawPicture();
}