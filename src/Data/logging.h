#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <string>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "settings.h"
#include "Data/dataExport.h"

void setLogFile(std::string simulationName, std::string dataPath);

#endif