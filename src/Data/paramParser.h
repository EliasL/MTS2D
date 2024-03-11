#ifndef PARAMPARSER_H
#define PARAMPARSER_H
#pragma once

#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <filesystem>
#include <stdexcept>
#include "alglibmisc.h"

// I am sorry about writing a custom parser.
// I tried to add a Yaml parser, but they were problamatic

struct Config
{
    // Simulation settings
    std::string name;
    int rows;
    int cols;
    int nrThreads;
    int seed;
    double plasticityEventThreshold;
    int usingPBC;
    // Loading settings
    double startLoad;
    double loadIncrement;
    double maxLoad;
    double noise;
    // Alglib settings
    int nrCorrections;
    double scale;
    double epsg;
    double epsf;
    double epsx;
    alglib::ae_int_t maxIterations;
    // Logging settings
    int showProgress;

    friend std::ostream &operator<<(std::ostream &os, const Config &config);
    std::string str() const;
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config getConf(std::string configFile);

#endif