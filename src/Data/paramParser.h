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

struct Config
{
    std::string name;
    int nx;
    int ny;
    int nrThreads;
    int seed;
    double startLoad;
    double loadIncrement;
    double maxLoad;
    double noise;
    int nrCorrections;
    double epsg;
    double epsf;
    double epsx;
    alglib::ae_int_t maxIterations;

    friend std::ostream &operator<<(std::ostream &os, const Config &config);
    std::string str() const;
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config getConf(std::string configFile);

#endif