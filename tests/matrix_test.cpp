#include "run/doctest.h"
#include "../src/matrix.h"  // Include the header for your Matrix class

TEST_CASE("Matrix Equality") {
    Matrix<int> matrix1(2, 2);
    Matrix<int> matrix2(2, 2);
    matrix1.fill(1);
    matrix2.fill(1);

    // Test if two equal matrices are indeed equal
    REQUIRE(matrix1 == matrix2);

    // Modify one element to make them unequal
    matrix1[1][1] = 2;
    REQUIRE(!(matrix1 == matrix2));
}

TEST_CASE("Matrix Addition") {
    Matrix<int> matrix1(2, 2);
    Matrix<int> matrix2(2, 2);
    matrix1.fill(1);
    matrix2.fill(2);

    // Perform matrix addition
    Matrix<int> result = matrix1 + matrix2;

    Matrix<int> expected(2, 2);
    expected.fill(3);

    // Test if the result is as expected
    REQUIRE(result == expected);
}

TEST_CASE("Matrix Multiplication") {
    Matrix<int> A(2, 3);
    Matrix<int> B(3, 2);

    // Initialize matrices A and B with values
    A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
    A[1][0] = 4; A[1][1] = 5; A[1][2] = 6;

    B[0][0] = 7; B[0][1] = 8;
    B[1][0] = 9; B[1][1] = 10;
    B[2][0] = 11; B[2][1] = 12;

    Matrix<int> result = A * B;

    Matrix<int> expected(2, 2);
    expected[0][0] = 58; expected[0][1] = 64;
    expected[1][0] = 139; expected[1][1] = 154;

    REQUIRE(result == expected);
}

TEST_CASE("Matrix Determinant") {
    // Test case 1: Determinant of a 1x1 matrix
    Matrix<int> matrix1x1(1, 1);
    matrix1x1[0][0] = 5;
    int det1x1 = matrix1x1.det();
    REQUIRE(det1x1 == 5);

    // Test case 2: Determinant of a 2x2 matrix
    Matrix<int> matrix2x2(2, 2);
    matrix2x2[0][0] = 1;
    matrix2x2[0][1] = 2;
    matrix2x2[1][0] = 3;
    matrix2x2[1][1] = 4;
    int det2x2 = matrix2x2.det();
    REQUIRE(det2x2 == -2);

    // Test case 3: Determinant of a 3x3 matrix
    Matrix<int> matrix3x3(3, 3);
    matrix3x3[0][0] = 2;
    matrix3x3[0][1] = 0;
    matrix3x3[0][2] = 1;
    matrix3x3[1][0] = 0;
    matrix3x3[1][1] = 1;
    matrix3x3[1][2] = 0;
    matrix3x3[2][0] = 1;
    matrix3x3[2][1] = 0;
    matrix3x3[2][2] = 2;
    int det3x3 = matrix3x3.det();
    REQUIRE(det3x3 == 3);
}
