#include "run/doctest.h"
#include "../src/matrix2x2.h"  // Include the header for your Matrix2x2 class

TEST_CASE("Matrix2x2 Initialization") {
    // Create a 2x2 matrix with known values for testing
    Matrix2x2<int> mat(1, 2, 3, 4);

    // Verify that the matrix elements are correct
    REQUIRE(mat[0][0] == 1);
    REQUIRE(mat[0][1] == 2);
    REQUIRE(mat[1][0] == 3);
    REQUIRE(mat[1][1] == 4);
}

TEST_CASE("Matrix2x2 Multiplication") {
    Matrix2x2<int> mat1(1, 2, 3, 4);
    Matrix2x2<int> mat2(5, 6, 7, 8);
    Matrix2x2<int> result = mat1 * mat2;

    // Verify the result of matrix multiplication
    REQUIRE(result[0][0] == 19);
    REQUIRE(result[0][1] == 22);
    REQUIRE(result[1][0] == 43);
    REQUIRE(result[1][1] == 50);
}

TEST_CASE("Matrix2x2 Determinant") {
    Matrix2x2<int> mat(1, 2, 3, 4);

    // Calculate the determinant
    int det = mat.det();

    // Verify the determinant value
    REQUIRE(det == -2);
}

TEST_CASE("Matrix2x2 Addition") {
    Matrix2x2<int> mat1(1, 2, 3, 4);
    Matrix2x2<int> mat2(5, 6, 7, 8);
    Matrix2x2<int> result = mat1 + mat2;

    // Verify the result of matrix addition
    REQUIRE(result[0][0] == 6);
    REQUIRE(result[0][1] == 8);
    REQUIRE(result[1][0] == 10);
    REQUIRE(result[1][1] == 12);
}

TEST_CASE("Matrix2x2 Transpose") {
    Matrix2x2<int> mat(1, 2, 3, 4);
    Matrix2x2<int> transposed = mat.transpose();

    // Verify the transpose of the matrix
    REQUIRE(transposed[0][0] == 1);
    REQUIRE(transposed[0][1] == 3);
    REQUIRE(transposed[1][0] == 2);
    REQUIRE(transposed[1][1] == 4);
}

TEST_CASE("Matrix2x2 Inverse") {
    // Define a small epsilon value
    const double epsilon = 1e-15;
    
    Matrix2x2<double> mat(2.0, 1.0, 1.0, 3.0);
    Matrix2x2<double> inverse = mat.inverse();

    // Verify the inverse of the matrix
    REQUIRE(inverse[0][0] - 3/5. < epsilon);
    REQUIRE(inverse[0][1] - -1/5. < epsilon);
    REQUIRE(inverse[1][0] - -1/5. < epsilon);
    REQUIRE(inverse[1][1] - 2/5. < epsilon);
}
