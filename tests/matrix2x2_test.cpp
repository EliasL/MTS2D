#include "run/doctest.h"
#include "../src/Matrix/matrix2x2.h"  // Include the header for your Matrix2x2 class
#include "../src/settings.h"

TEST_CASE("Matrix2x2 Initialization") {
    // Create a 2x2 matrix with known values for testing
    Matrix2x2<int> mat(1, 2, 3, 4);

    // Verify that the matrix elements are correct
    REQUIRE(mat[0][0] == 1);
    REQUIRE(mat[0][1] == 2);
    REQUIRE(mat[1][0] == 3);
    REQUIRE(mat[1][1] == 4);
}


TEST_CASE("Matrix2x2 swapCols Test") {
    Matrix2x2<int> mat(1, 2, 3, 4);
    mat.swapCols();
    // TODO double check this
    // Are rows rows? Or are they columns? Make some documentation
    CHECK(mat[0][0] == 2);
    CHECK(mat[0][1] == 1);
    CHECK(mat[1][0] == 4);
    CHECK(mat[1][1] == 3);
}

TEST_CASE("Matrix2x2 flip Test") {
    Matrix2x2<int> mat(1, -2, 3, 4);
    mat.flip(0, 0);
    mat.flip(0, 1);
    mat.flip(1, 0);
    mat.flip(1, 1);
    
    CHECK(mat[0][0] == -1);
    CHECK(mat[0][1] == 2);
    CHECK(mat[1][0] == -3);
    CHECK(mat[1][1] == -4);
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

TEST_CASE("Matrix2x2 lag_m1 Test") {
    Matrix2x2<int> mat(1, 2, 3, 4);
    mat.lag_m1();
    
    CHECK(mat[0][0] == 1);
    CHECK(mat[0][1] == 2);
    CHECK(mat[1][0] == -3);
    CHECK(mat[1][1] == -4);
}

TEST_CASE("Matrix2x2 lag_m2 Test") {
    Matrix2x2<int> mat(1, 2, 3, 4);
    mat.lag_m2();
    
    CHECK(mat[0][0] == 2);
    CHECK(mat[0][1] == 1);
    CHECK(mat[1][0] == 4);
    CHECK(mat[1][1] == 3);
}

TEST_CASE("Matrix2x2 lag_m3 Test") {
    Matrix2x2<int> mat(1, 2, 3, 4);
    mat.lag_m3();
    
    CHECK(mat[0][0] == 3);
    CHECK(mat[0][1] == 1);
    CHECK(mat[1][0] == 7);
    CHECK(mat[1][1] == 1);
}

TEST_CASE("Test conjugate function") {
    Matrix2x2<double> a;
    a[0][0] = 1;
    a[0][1] = 2;
    a[1][0] = 3;
    a[1][1] = 4;

    Matrix2x2<double> g;
    g[0][0] = 0.5;
    g[0][1] = 0.2;
    g[1][0] = -0.3;
    g[1][1] = 0.8;

    Matrix2x2<double> result = a.conjugate(g);
    //https://www.wolframalpha.com/input?i=%7B%7B0.5%2C+0.2%7D%2C+%7B-0.3%2C+0.8%7D%7D++%7B%7B1%2C2%7D%2C%7B3%2C4%7D%7D++inverse+%7B%7B0.5%2C+0.2%7D%2C+%7B-0.3%2C+0.8%7D%7D+
    
    double expected00 = 3.08696;
    double expected01 = 1.47826;
    double expected10 = 5.34783 ;
    double expected11 = 1.91304;

    CHECK(result[0][0] == doctest::Approx(expected00));
    CHECK(result[0][1] == doctest::Approx(expected01));
    CHECK(result[1][0] == doctest::Approx(expected10));
    CHECK(result[1][1] == doctest::Approx(expected11));
}

TEST_CASE("Test orth_conjugate function") {
    Matrix2x2<double> a;
    a[0][0] = 1;
    a[0][1] = 2;
    a[1][0] = 3;
    a[1][1] = 4;

    Matrix2x2<double> g = Matrix2x2<double>::rotation_matrix(3.0/4.0*M_PI);

    Matrix2x2<double> result = a.orth_conjugate(g);
    //https://www.wolframalpha.com/input?i=%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D++%7B%7B1%2C2%7D%2C%7B3%2C4%7D%7D++transpose+%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D+where+a%3D3pi%2F4+    double expected00 = 5;
    double expected00 = 5;
    double expected01 = 1;
    double expected10 = 2 ;
    double expected11 = 0;

    CHECK(result[0][0] == doctest::Approx(expected00));
    CHECK(result[0][1] == doctest::Approx(expected01));
    CHECK(result[1][0] == doctest::Approx(expected10));
    CHECK(result[1][1] == doctest::Approx(expected11));
}

TEST_CASE("Test isOrthogonal function") {
    // Create an orthogonal matrix and test 'isOrthogonal' function.
    Matrix2x2<double> g;
    g[0][0] = 1;
    g[0][1] = 0;
    g[1][0] = 0;
    g[1][1] = -1;

    CHECK(g.isOrthogonal() == true);

    g = Matrix2x2<double>(3);
    CHECK(g.isOrthogonal() == false);

    // 3 is a random theta. No special meaning.
    g = Matrix2x2<double>::reflection_matrix(2.0);
    CHECK(g.isOrthogonal() == true);
}

TEST_CASE("Test sym_orth_conjugate function") {
    // Create a symmetric matrix and an orthogonal matrix and test 'sym_orth_conjugate' function.
    // Make sure to define your expected result based on the known computation of the symmetric orthogonal conjugate.
    // Replace the following values with the expected values.
    Matrix2x2<double> a;
    a[0][0] = 1;
    a[0][1] = 2;
    a[1][0] = 2; // Symmetric element
    a[1][1] = 3;

    Matrix2x2<double> g = Matrix2x2<double>::rotation_matrix(3.0/4.0*M_PI);
    Matrix2x2<double> result = a.sym_orth_conjugate(g);

    //https://www.wolframalpha.com/input?i=%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D++%7B%7B1%2C2%7D%2C%7B2%2C3%7D%7D++inverse+%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D+where+a%3D3pi%2F4+
    double expected00 = 4;
    double expected01 = 1;
    double expected10 = 1;
    double expected11 = 0;

    CHECK(result[0][0] == doctest::Approx(expected00));
    CHECK(result[0][1] == doctest::Approx(expected01));
    CHECK(result[1][0] == doctest::Approx(expected10));
    CHECK(result[1][1] == doctest::Approx(expected11));
}

