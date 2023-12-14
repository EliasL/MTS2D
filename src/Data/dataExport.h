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

#include "../Mesh/mesh.h"
#include "../Matrix/matrix.h"

bool create_directory_if_not_exists(const std::string &path)
{
    std::vector<std::string> dirs;
    std::stringstream ss(path);
    std::string item;
    while (getline(ss, item, '/'))
    {
        if (!item.empty())
        {
            dirs.push_back(item);
        }
    }

    std::string current_path;
    for (const auto &dir : dirs)
    {
        current_path += dir + "/";
#if defined(_WIN32)
        // Windows does not have a built-in function for recursive directory creation
        // Here you might want to implement a loop that creates each directory in the path
        if (_mkdir(current_path.c_str()) != 0 && errno != EEXIST)
        {
            std::cerr << "Error creating directory '" << current_path << "': " << strerror(errno) << std::endl;
            return false;
        }
#else
        struct stat info;
        if (stat(current_path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR))
        {
            if (mkdir(current_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0 && errno != EEXIST)
            {
                std::cerr << "Error creating directory '" << current_path << "': " << strerror(errno) << std::endl;
                return false;
            }
        }
#endif
    }

    return true;
}

std::string getFilePath(std::string fileName,
                        std::string fileType = "vtu",
                        std::string folder = "output/")
{
    std::string fileNameWithType = fileName + '.' + fileType;
    std::string directory = folder + fileType + '/';
        // Ensure the directory exists
    if (!create_directory_if_not_exists(directory))
    {
        std::cerr << "Failed to create directory: " << directory << std::endl;
    }

    return directory + fileNameWithType;
}

std::ofstream getWriteFile(std::string fileName,
                           std::string fileType = "vtu",
                           std::string folder = "output/")
{
    std::string directory = getFilePath(fileName, folder, fileType);

    std::ofstream filestr(directory + fileName); // Open the file with the specified mode
    if (!filestr.is_open())
    {
        std::cerr << "Failed to open file for writing: " << fileName << std::endl;
    }
    return filestr;
}

void write_to_vtu(Mesh &mesh, const std::string &fileName = "data")
{
    const int dim = 3;
    const int cell_size = 3;
    int nrNodes = mesh.nodes.data.size();
    int nrElements = mesh.nrElements;
    std::string filePath = getFilePath(fileName);

    std::vector<double> points(nrNodes * dim, 0);
    std::vector<double> force(nrNodes * dim, 0);
    std::vector<int> elements(nrElements * cell_size);
    std::vector<double> energy(nrElements);

    leanvtk::VTUWriter writer;

    // Transfer data
    for (int i = 0; i < nrNodes; ++i)
    {
        points[i * 3 + 0] = mesh.nodes.data[i].x;
        points[i * 3 + 1] = mesh.nodes.data[i].y;
        points[i * 3 + 2] = 0;
        force[i * 3 + 0] = mesh.nodes.data[i].f_x;
        force[i * 3 + 1] = mesh.nodes.data[i].f_y;
        force[i * 3 + 2] = 0;
    }

    for (int i = 0; i < nrElements; ++i)
    {
        elements[i * 3 + 0] = mesh.elements[i].n1->id.i;
        elements[i * 3 + 1] = mesh.elements[i].n2->id.i;
        elements[i * 3 + 2] = mesh.elements[i].n3->id.i;
        energy[i] = mesh.elements[i].energy;
    }

    // connect data to writer
    writer.add_cell_scalar_field("energy_field", energy);
    writer.add_vector_field("stress_field", force, dim);

    // write data
    writer.write_surface_mesh(filePath, dim, cell_size, points, elements);
}

void write_to_legacy_vtk(Mesh &mesh, const std::string &fileName = "data")
{
    int nrNodes = mesh.nodes.data.size();
    int nrElements = mesh.nrElements;
    std::ofstream filestr = getWriteFile(fileName, "vtk");

    // Write VTK header
    filestr << "# vtk DataFile Version 3.0\n";
    filestr << "Unstructured Grid Example\n";
    filestr << "ASCII\n";
    filestr << "DATASET UNSTRUCTURED_GRID\n";

    filestr << "POINTS " << nrNodes << " float\n";
    // Write point data
    for (int i = 0; i < nrNodes; ++i)
    {
        float x = static_cast<float>(mesh.nodes.data[i].x);
        float y = static_cast<float>(mesh.nodes.data[i].y);
        float z = 0.0f; // Assuming 2D data, z-coordinate is 0
        filestr << x << " " << y << " " << z << "\n";
    }

    filestr << "CELLS " << nrElements << " " << 4 * nrElements << "\n";
    for (int i = 0; i < nrElements; ++i)
    {
        int n1Id = mesh.elements[i].n1->id.i;
        int n2Id = mesh.elements[i].n2->id.i;
        int n3Id = mesh.elements[i].n3->id.i;
        filestr << "3 " << n1Id << " " << n2Id << " " << n3Id << "\n";
    }

    filestr << "CELL_TYPES " << nrElements << "\n";
    for (int i = 0; i < nrElements; ++i)
    {
        filestr << "5\n"; // 5 is the VTK cell type for a triangle
    }

    // Write point data (scalars, vectors, etc.)
    filestr << "POINT_DATA " << nrNodes << "\n";
    filestr << "SCALARS energy float 1\n";
    filestr << "LOOKUP_TABLE default\n";

    for (int i = 0; i < nrNodes; ++i)
    {
        float scalarValue = static_cast<float>(mesh.elements[i].energy); // Example scalar
        filestr << scalarValue << "\n";
    }

    filestr.close();
}

void write_to_xyz(Mesh &mesh, std::string file_name = "data")
{

    int n = mesh.nodes.data.size();
    std::string directory = "output/ovito/";
    std::string filename = directory + file_name + ".xyz";

    // Ensure the directory exists
    if (!create_directory_if_not_exists(directory))
    {
        std::cerr << "Failed to create directory: " << directory << std::endl;
        return;
    }

    std::ofstream filestr(filename.c_str());
    if (!filestr.is_open())
    {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    filestr << n << "\n";
    filestr << " "
            << "\n";

    for (int i = 0; i < n; ++i)
    {

        // C is the metrics of the cell
        double sqroot_detC = sqrt(mesh.elements[i].C.det());

        filestr << std::scientific << std::setprecision(16)
                << mesh.nodes.data[i].x << " "
                << mesh.nodes.data[i].y << " "
                << mesh.nodes.data[i].f_x << " "
                << mesh.nodes.data[i].f_y << " "
                << mesh.elements[i].energy << " "
                << mesh.elements[i].C[1][0] << " "
                << mesh.elements[i].C[0][0] << " "
                << mesh.elements[i].C[1][1] << " "
                << mesh.elements[i].P[1][0] << " "
                << mesh.elements[i].P[0][0] << " "
                << mesh.elements[i].P[1][1] << " "
                //<< contraction(c.stress[i][k],setnew) << " "

                << sqroot_detC << " "
                << mesh.load
                << std::endl;
    }

    filestr.close();
}

#endif