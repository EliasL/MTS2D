#include "logging.h"

void setLogFile(std::string simulationName)
{

    std::string path = getOutputPath(simulationName);
    std::string filePath = path + simulationName + ".log";
    // We set the logging settings
    auto file_logger = spdlog::basic_logger_mt(LOGNAME, filePath);
    spdlog::set_default_logger(file_logger);
    spdlog::info("Starting simulation.");
}