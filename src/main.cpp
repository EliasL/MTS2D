#include <iostream>

#include "spdlog/spdlog.h"

#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

int main(){ 
    
    // We fix the random seed to get reproducable results
    srand(0);
    run_simulation();
    // drawPicture();
}