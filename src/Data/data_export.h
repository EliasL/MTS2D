#ifndef data_export_H
#define data_export_H
#include "Data/param_parser.h"
#pragma once

#include <fstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "../Mesh/mesh.h"

// Get date
std::string getCurrentDate();

// If no outputPath is provided, we try to automatically search for a existing
// path
std::string findOutputPath();
std::string searchForConfig(std::string dumpPath);
/*
Each simulation run should take place in its own folder. The folder will have
two subfolders containing the raw data and frames. A smaller cvs file will
contain a small amount of processed data, one line for each frame, as opposed
to one file per frame as done inside the data folder.
*/
std::string getOutputPath(const std::string &name, const std::string &dataPath);
std::string getDataPath(const std::string &name, const std::string &dataPath);
std::string getMinDataSubFolder(const Mesh mesh);
std::string getFramePath(const std::string &name, const std::string &dataPath);
std::string getDumpPath(const std::string &name, const std::string &dataPath);
std::string getBackupPath(const std::string &name, const std::string &dataPath);

// Generates a name based on the settings provided
std::string makeFileName(const Mesh &mesh, std::string name,
                         std::string dataPath);

// Clears a subfolder. It only clears .vtu and .pvd files for safety.
// If you want to delete the entire outputfolder, do it manually.
void clearOutputFolder(std::string name, std::string dataPath);

// Each frame (load step) can be saved to a seperate Vtu file
std::string writeMeshToVtu(const Mesh &mesh, std::string folderName,
                           std::string dataPath, std::string fileName = "",
                           bool minimizationStep = false);

// Duplicated the config file into the output
void saveConfigFile(std::string configFile, std::string dataPath);
void saveConfigFile(Config conf, std::string dataPath);

void writeLineToCsv(std::ofstream &file,
                    const std::vector<std::string> &strings);
void writeLineToCsv(std::ofstream &file, const std::vector<double> &values);
std::vector<std::string> getCsvCols();

// Forward declaration of simulation class
class Simulation;

// The averaged values of each frame can be saved to a single cvs file
// The first row of the cvs file should indicate the name of the columns
// eg. Frame nr, Avg. energy, Avg. Stress, Nr. dislocations
// This function is also used to create other
std::ofstream initCsvFile(const std::string &folderName,
                          const std::string &dataPath, const Simulation &s,
                          const std::string subFolder = "");
// When a simulation is resumed from a dump, unless the program was stopped
// right after the dump was created, the csv file will have lines that need to
// be overwritten
void trimCsvFile(const std::string &file, const Simulation &s);
std::vector<std::string> getStringVector(const Simulation &s);
void writeToCsv(std::ofstream &file, const Simulation &s);
void writeCsvCols(std::ofstream &file);
// returns true if header was written
bool insertHeaderIfNeeded(const std::string &filename);

/**
 * Finds all files of specified type and creates a .pvd collection
 * The files must be on the form name.N.x, where N is an integer number,
 * and x is the extension of the file.
 * const std::string folderPath         The path to a folder with the files
 *                                      that the collection should link
 * const std::string extension          The extension of the files to be linked
 * const std::vector<double>& timestep  An optional parameter to set the
 * timestep of each frame.
 */

void createCollection(
    const std::string &folderPath, const std::string &destination,
    const std::string &regexPattern = "", const std::string &extension = ".vtu",
    const std::vector<double> &timestep = std::vector<double>());

// Reads the last line of the cvs file and checks whether or not maxLoad is
// reached or not
bool simulationAlreadyComplete(const std::string &name,
                               const std::string &dataPath, double maxLoad);

#endif