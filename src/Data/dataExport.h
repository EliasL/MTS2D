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
#include "lean_vtk.h"

#if defined(_WIN32)
#include <direct.h> // For _mkdir on Windows
#else
#include <unistd.h>
#endif

#include "settings.h"
#include "../Mesh/mesh.h"
#include "../Matrix/matrix.h"

// Generates a name based on the settings provided
std::string makeFileName(const Mesh &mesh);

// Creates a folder inside the output folder, and creates two folders
// inside for data and frames
void createDataFolder();

// Clears a subfolder. It only clears .vtu and .pvd files for safety.
// If you want to delete the entire outputfolder, do it manually.
void clearOutputFolder();

// Sets the config of the logger to output the log file into the desired folder
void setLoggingOutput();

void writeToVtu(Mesh &mesh, bool automaticNumbering = true);

#endif