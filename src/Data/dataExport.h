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


#if defined(_WIN32)
#include <direct.h> // For _mkdir on Windows
#else
#include <unistd.h>
#endif

#include "../Surface/surface.h"
#include "../Matrix/matrix.h"

bool create_directory_if_not_exists(const std::string &path) {
    std::vector<std::string> dirs;
    std::stringstream ss(path);
    std::string item;
    while (getline(ss, item, '/')) {
        if (!item.empty()) {
            dirs.push_back(item);
        }
    }

    std::string current_path;
    for (const auto &dir : dirs) {
        current_path += dir + "/";
        #if defined(_WIN32)
        // Windows does not have a built-in function for recursive directory creation
        // Here you might want to implement a loop that creates each directory in the path
        if (_mkdir(current_path.c_str()) != 0 && errno != EEXIST) {
            std::cerr << "Error creating directory '" << current_path << "': " << strerror(errno) << std::endl;
            return false;
        }
        #else
        struct stat info;
        if (stat(current_path.c_str(), &info) != 0 || !(info.st_mode  &S_IFDIR)) {
            if (mkdir(current_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0 && errno != EEXIST) {
                std::cerr << "Error creating directory '" << current_path << "': " << strerror(errno) << std::endl;
                return false;
            }
        }
        #endif
    }

    return true;
}


// void write_to_vtk(const struct conf_stru &c, int t){

//         int id;
//         string filename;

//         int nx =  c.nx;
//         int ny =  c.ny;
//         int n =   c.n;
//         int nt=   c.myVector.size();

//         filename = get_outputPaht() + "/dir_vtk/vtk_" + IntToStr(t) + ".vtk";
//         std::ofstream par;
//         par.open(filename.c_str()); // nom du fichier qui contient le maillage
//                                     // 
//         par <<"# vtk DataFile Version 1.0" << std::endl;
//         par <<"2D Unstructured Mesh of Linear Triangles" << std::endl;
//         par <<"ASCII" << std::endl;
//         par <<" " << std::endl;

//         par <<"DATASET UNSTRUCTURED_GRID" << std::endl;
//         //// representation du solide deformé
//         par <<"POINTS " << n<< " float" << std::endl;
//         for (int id=0;id<n;id++) par<< c.p[id].x <<" "<< c.p[id].y <<"  "<< 0 <<std::endl;
//         par <<"CELLS " <<  nt<< " " << (4)*nt<<std::endl; // 4 car 3 dim + 1
//         for (const auto &element : c.myVector) {
//                 int i = element.first;
//                 int k = element.second;
//                 struct punto_stru p = c.p[i];
//                 par<<"3 "<< i <<" "<< p.nn[k]<<" "<<p.nn[k+1]<<std::endl;
//         };
//         par <<" " <<std::endl;
//         par <<"CELL_TYPES " << nt<<std::endl;
//         for (int id=0;id<nt;id++) par<<"5"<<std::endl;
//         /////////////////////Visualisation des champs sur la configuration deformée
//         par <<" " <<std::endl;
//         par <<"POINT_DATA " << n <<std::endl;
//         /// energie potentielle
//         par << "SCALARS Energy float" << std::endl;
//         par << "LOOKUP_TABLE default " << std::endl;
         
//         for (int id=0;id<n;id++) par << c.energy_local[id][0] << std::endl;

//         par <<" " <<std::endl;
//         par <<"CELL_DATA " << nt <<std::endl;
//         par <<"TENSORS Stress float " <<std::endl;
// //      par << "LOOKUP_TABLE default " << std::endl;
//         for (const auto &element : c.myVector) {
//                 int i = element.first;
//                 int k = element.second;
         
//                 par << c.stress[i][k].c11 << " " << c.stress[i][k].c12 << " " << 0 << std::endl;
//                 par << c.stress[i][k].c12 << " " << c.stress[i][k].c22 << " " << 0 << std::endl;
//                 par << 0 << " " << 0 << " " << 0 << std::endl;
//                 par << " " << std::endl;
//         };

//         //// deviateur
// //      par <<"SCALARS Deviator float"<<std::endl;
// //      par <<"LOOKUP_TABLE default " <<std::endl;
// //      for (int id=0;id<n;id++) {par<<  devVMdes[][id]<<std::endl; }
// };

void write_to_a_ovito_file(Mesh &mesh, std::string file_name = "data"){

    int n = mesh.nodes.data.size();
    std::string directory = "output/ovito/";
    std::string filename = directory + file_name + ".xyz";

    // Ensure the directory exists
    if (!create_directory_if_not_exists(directory)) {
        std::cerr << "Failed to create directory: " << directory << std::endl;
        return;
    }

    std::ofstream filestr(filename.c_str());
    if (!filestr.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    filestr << n << "\n";
    filestr << " " << "\n";

	for (int i = 0; i<n; ++i){

        // C is the metrics of the cell
		double sqroot_detC = sqrt(mesh.cells[i].C.det());

		filestr << std::scientific << std::setprecision(16)
			<< mesh.nodes.data[i].x << " "
			<< mesh.nodes.data[i].y << " "
			<< mesh.nodes.data[i].f_x << " "
			<< mesh.nodes.data[i].f_y << " "
			<< mesh.cells[i].energy  << " "
			<< mesh.cells[i].C[1][0] << " "
			<< mesh.cells[i].C[0][0] << " "
			<< mesh.cells[i].C[1][1] << " "
			<< mesh.cells[i].P[1][0] << " "
			<< mesh.cells[i].P[0][0] << " "
			<< mesh.cells[i].P[1][1] << " " 
			//<< contraction(c.stress[i][k],setnew) << " "

			<< sqroot_detC << " "
			<< mesh.load
			<< std::endl;
	}

 	filestr.close();
}

#endif