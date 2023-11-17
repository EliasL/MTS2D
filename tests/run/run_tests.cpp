#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

// Log settings
//#define ELPP_DISABLE_LOGS
//#define ELPP_DISABLE_DEBUG_LOGS
#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

void logging_config(){
    // Load configuration from file
    el::Configurations conf("../libs/easylogging++/logging.conf");
    // Reconfigure single logger
    el::Loggers::reconfigureLogger("default", conf);
    // Actually reconfigure all loggers instead
    el::Loggers::reconfigureAllLoggers(conf);
    // Now all the loggers will use configuration from file
}


int main(int argc, char* argv[]) {
    logging_config();

    doctest::Context context;
    context.applyCommandLine(argc, argv);
    int res = context.run(); // Run tests
    return res;
}