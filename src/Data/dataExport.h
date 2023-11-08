#ifndef DATAEXPORT_H
#define DATAEXPORT_H
#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "../Grid/grid2D.h"
#include "../Matrix/matrix.h"


// void write_to_vtk(const struct conf_stru& c, int t){

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
//         par <<"2D Unstructured Grid of Linear Triangles" << std::endl;
//         par <<"ASCII" << std::endl;
//         par <<" " << std::endl;

//         par <<"DATASET UNSTRUCTURED_GRID" << std::endl;
//         //// representation du solide deformé
//         par <<"POINTS " << n<< " float" << std::endl;
//         for (int id=0;id<n;id++) par<< c.p[id].x <<" "<< c.p[id].y <<"  "<< 0 <<std::endl;
//         par <<"CELLS " <<  nt<< " " << (4)*nt<<std::endl; // 4 car 3 dim + 1
//         for (const auto& element : c.myVector) {
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
//         for (const auto& element : c.myVector) {
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

void write_to_a_ovito_file(Grid& g, std::string file_name = "data"){

	int n =   g.nodes.data.size();

	std::string filename;

	filename = "output/ovito/" + file_name + ".xyz";

	std::ofstream filestr;
	filestr.open(filename.c_str());

	filestr << n << std::endl;
	filestr << " " << std::endl;

	for (int i = 0; i<n; ++i){

        // C is the metrics of the cell
		double sqroot_detC = sqrt(g.cells[i].C.det());

		filestr << std::scientific << std::setprecision(16)
			<< g.nodes.data[i].x << " "
			<< g.nodes.data[i].y << " "
			<< g.cells[i].energy  << " "
			<< g.cells[i].C[1][0] << " "
			<< g.cells[i].C[0][0] << " "
			<< g.cells[i].C[1][1] << " "
			<< g.cells[i].P[1][0] << " "
			<< g.cells[i].P[0][0] << " "
			<< g.cells[i].P[1][1] 
			//<< contraction(c.stress[i][k],setnew) << " "

			<< sqroot_detC << " "
			<< g.load
			<< std::endl;
	}

 	filestr.close();

//****************************************************************************************



}

#endif