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
#include <boost/serialization/access.hpp>

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
    std::string scenario;
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
    int maxIterations;
    // Logging settings
    int showProgress;

    friend std::ostream &operator<<(std::ostream &os, const Config &config);
    std::string str() const;

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version);
};

std::map<std::string, std::string> parseParams(const std::string &filename);

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap);

Config getConf(std::string configFile);

#endif

template <class Archive>
void Config::serialize(Archive &ar, const unsigned int version)
{
    ar & name;
    ar & rows;
    ar & cols;
    ar & nrThreads;
    ar & seed;
    ar & plasticityEventThreshold;
    ar & scenario;
    ar & startLoad;
    ar & loadIncrement;
    ar & maxLoad;
    ar & noise;
    ar & nrCorrections;
    ar & scale;
    ar & epsg;
    ar & epsf;
    ar & epsx;
    ar & maxIterations;
    ar & showProgress;
}