#include "paramParser.h"

std::string Config::str() const
{
    std::stringstream ss;
    ss << *this;
    return ss.str();
}

std::ostream &operator<<(std::ostream &os, const Config &config)
{
    os << "Name: " << config.name << "\n"
       << "rows, cols: " << config.rows << ", " << config.cols << "\n"
       << "boundary conditions: " << (config.usingPBC ? "PBC" : "NPBC") << "\n"
       << "scenario: " << config.scenario << "\n"
       << "nrThreads: " << config.nrThreads << "\n"
       << "seed: " << config.seed << "\n"
       << "plasticityEventThreshold: " << config.plasticityEventThreshold << "\n"
       << "startLoad, loadIncrement, maxLoad: "
       << config.startLoad << ", " << config.loadIncrement << ", " << config.maxLoad << "\n"
       << "noise: " << config.noise << "\n"
       << "nrCorrections: " << config.nrCorrections << "\n"
       << "scale: " << config.scale << "\n"
       << "epsg, epsf, epsx: "
       << config.epsg << ", " << config.epsf << ", " << config.epsx << "\n"
       << "maxIterations: " << config.maxIterations << "\n"
       << "showProgress: " << config.showProgress << "\n";
    return os;
}

namespace fs = std::filesystem;

std::map<std::string, std::string> parseParams(const std::string &filename)
{
    std::map<std::string, std::string> config;
    std::ifstream file(filename);
    std::string line;

    // Extract the file name from the path and assign it to "name"
    fs::path filePath(filename);
    if (!file)
    {                                                            // Check if the file was successfully opened
        throw std::runtime_error("File not found: " + filename); // Throw an error if not
    }
    config["name"] = filePath.stem().string();

    while (std::getline(file, line))
    {
        // Remove comments (anything after '#')
        std::string::size_type commentPos = line.find('#');
        if (commentPos != std::string::npos)
        {
            line = line.substr(0, commentPos);
        }

        // Trim the line for whitespace before parsing
        line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());

        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '='))
        {
            std::string value;
            if (std::getline(is_line, value))
            {
                config[key] = value;
            }
        }
    }

    config["name"] = filePath.stem().string();

    return config;
}

// Function to initialize Config from a map
Config initializeConfig(const std::map<std::string, std::string> &configMap)
{
    Config config;

    // Initialize other variables directly from configMap with the appropriate conversions

    config.name = configMap.at("name");
    config.rows = std::stoi(configMap.at("rows"));
    config.cols = std::stoi(configMap.at("cols"));
    config.usingPBC = std::stoi(configMap.at("usingPBC")); // Note that this is a bool, not int
    config.scenario = configMap.at("scenario");
    config.nrThreads = std::stoi(configMap.at("nrThreads"));
    config.seed = std::stoi(configMap.at("seed"));
    config.plasticityEventThreshold = std::stod(configMap.at("plasticityEventThreshold"));
    config.startLoad = std::stod(configMap.at("startLoad"));
    config.loadIncrement = std::stod(configMap.at("loadIncrement"));
    config.maxLoad = std::stod(configMap.at("maxLoad"));
    config.noise = std::stod(configMap.at("noise"));
    config.nrCorrections = std::stoi(configMap.at("nrCorrections"));
    config.scale = std::stod(configMap.at("scale"));
    config.epsg = std::stod(configMap.at("epsg"));
    config.epsf = std::stod(configMap.at("epsf"));
    config.epsx = std::stod(configMap.at("epsx"));
    config.maxIterations = std::stoi(configMap.at("maxIterations"));
    config.showProgress = std::stoi(configMap.at("showProgress"));

    return config;
}

Config getConf(std::string configFile)
{
    auto confMap = parseParams(configFile);
    return initializeConfig(confMap);
}
