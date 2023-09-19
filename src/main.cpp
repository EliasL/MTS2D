#include "matrix.h"
#include "matrix2x2.h"
#include "grid2D.h"
#include "settings.h"
#include <iostream>

int main() {
    Matrix<int> intMatrix(3, 3);
    Matrix<double> doubleMatrix(3, 3);

    doubleMatrix[2][1]=1.1;
    doubleMatrix[0][2]=3;
    doubleMatrix[1][2]=6;
    doubleMatrix[2][0]=7;
    doubleMatrix[2][0]=7;
    doubleMatrix[0][0]=1;
    doubleMatrix[2][2]=2;

    std::cout << doubleMatrix << std::endl;
    // Use intMatrix and doubleMatrix with their respective types

    // Create a grid
    Grid g(4, 4);
    int x = 1;

    // Change the value of the atom at row 3, column 2
    g.atoms[3][2].x = 5;

    Matrix2x2<int> mat(0);
    mat[0][0] = 1;
    mat[0][1] = 2;
    mat[1][0] = 3;
    mat[1][1] = 4;
    mat[1][1] = 4;

}