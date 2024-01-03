#include <iostream>


// Log settings
//#define ELPP_DISABLE_LOGS
//#define ELPP_DISABLE_DEBUG_LOGS
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "Simulation/simulation.h"
#include "Data/dataExport.h"

void logging_config(){
    // Load configuration from file
    el::Configurations conf("../libs/easylogging++/logging.conf");
    // Reconfigure single logger
    el::Loggers::reconfigureLogger("default", conf);
    // Actually reconfigure all loggers instead
    el::Loggers::reconfigureAllLoggers(conf);
    // Now all the loggers will use configuration from file
}

int main(){ 
    clearOutputFolder();
    logging_config();
    run_simulation();
}