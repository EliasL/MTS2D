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
#include "lean_vtk.h"

#if defined(_WIN32)
#include <direct.h> // For _mkdir on Windows
#else
#include <unistd.h>
#endif

#include "settings.h"
#include "../Mesh/mesh.h"
#include "../Matrix/matrix.h"

// Clears a subfolder. It only clears .vtu and .pvd files for safety.
// If you want to delete the entire outputfolder, do it manually.
void clearOutputFolder(const std::string& folder = DEFAULTSUBFOLDER);

void writeToVtu(Mesh &mesh, std::string fileName = "data",
                  bool automaticNumbering = true);

void write_to_legacy_vtk(Mesh &mesh, const std::string &fileName = "data");

void write_to_xyz(Mesh &mesh, std::string file_name = "data");

#endif