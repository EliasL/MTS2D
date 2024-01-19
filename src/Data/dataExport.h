#ifndef DATAEXPORT_H
#define DATAEXPORT_H
#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <sys/stat.h>
#include <iterator>
#include <algorithm>
#include <cstring>
#include <vector>
#include <chrono>
#include <unistd.h>
#include "lean_vtk.h"

#include "settings.h"
#include "../Mesh/mesh.h"
#include "../Matrix/matrix.h"

/*
Each simulation run should take place in its own folder. The folder will have
two subfolders containing the raw data and frames. A smaller cvs file will
contain a small amount of processed data, one line for each frame, as opposed
to one file per frame as done inside the data folder.
*/

std::string getOutputPath(std::string name, std::string dataPath);
std::string getDataPath(std::string name, std::string dataPath);
std::string getFramePath(std::string name, std::string dataPath);

// Generates a name based on the settings provided
std::string makeFileName(const Mesh &mesh, std::string name, std::string dataPath);

// Creates a folder inside the output folder, and creates two folders
// inside for data and frames
void createDataFolder(std::string name, std::string dataPath);

// Clears a subfolder. It only clears .vtu and .pvd files for safety.
// If you want to delete the entire outputfolder, do it manually.
void clearOutputFolder(std::string name, std::string dataPath);

// Each frame (load step) can be saved to a seperate Vtu file
void writeToVtu(Mesh &mesh, std::string folderName, std::string dataPath,
                bool automaticNumbering = true);

// The averaged values of each frame can be saved to a single cvs file
// The first row of the cvs file should indicate the name of the columns
// eg. Frame nr, Avg. energy, Avg. Stress, Nr. dislocations
void writeLineToCsv(std::vector<std::string> &names, std::string folderName,
                    std::string dataPath);
void writeLineToCsv(std::vector<double> &values, std::string folderName,
                    std::string dataPath);

// Calculates the appropreate values and writes one line to the csv file
// Remember to use isFirstLine to write the names of the columns
void writeMeshToCsv(Mesh &mesh, std::string folderName, std::string dataPath,
                    bool isFirstLine = false);

#endif