#include "run/doctest.h"
#include "../src/Matrix/matrix2x2.h" // Include the header for your Matrix2x2 class
#include "../src/settings.h"

TEST_CASE("Matrix2x2 Initialization")
{
    // Create a 2x2 matrix with known values for testing
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});

    // Verify that the matrix elements are correct
    REQUIRE(mat[0][0] == 1);
    REQUIRE(mat[0][1] == 2);
    REQUIRE(mat[1][0] == 3);
    REQUIRE(mat[1][1] == 4);
}

TEST_CASE("Matrix2x2 swapCols Test")
{
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> ans({{2, 1},
                        {4, 3}});

    mat.swapCols();
    REQUIRE(mat == ans);
}

TEST_CASE("Matrix2x2 flip Test")
{
    Matrix2x2<int> mat({{1, -2},
                        {3, 4}});
    Matrix2x2<int> ans({{-1, 2},
                        {-3, -4}});
    mat.flip(0, 0);
    mat.flip(0, 1);
    mat.flip(1, 0);
    mat.flip(1, 1);

    REQUIRE(mat == ans);
}

TEST_CASE("Matrix2x2 set col")
{
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> ans({{5, 6},
                        {7, 8}});
    std::array<int, 2> c1 = {5, 7};
    std::array<int, 2> c2 = {6, 8};

    mat.setCol(c1, 0);
    mat.setCol(c2, 1);
    REQUIRE(mat == ans);
}

TEST_CASE("Matrix2x2 set cols")
{
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> ans({{5, 6},
                        {7, 8}});
    std::array<int, 2> c1 = {5, 7};
    std::array<int, 2> c2 = {6, 8};

    mat.setCols(c1, c2);
    REQUIRE(mat == ans);
}

TEST_CASE("Matrix2x2 Multiplication")
{
    Matrix2x2<int> mat1({{1, 2},
                         {3, 4}});
    Matrix2x2<int> mat2({{5, 6},
                         {7, 8}});
    Matrix2x2<int> result1 = mat1 * mat2;
    Matrix2x2<int> result2 = mat1 * mat2 * mat1;

    // https://www.wolframalpha.com/input?i=%7B%7B1%2C+2%7D%2C+%7B3%2C+4%7D%7D+%7B%7B5%2C+6%7D%2C+%7B7%2C+8%7D%7D
    //  Verify the result of matrix multiplication
    Matrix2x2<int> ans1({{19, 22},
                         {43, 50}});
    REQUIRE(result1 == ans1);

    // https://www.wolframalpha.com/input?i=+%7B%7B1%2C+2%7D%2C+%7B3%2C+4%7D%7D%7B%7B5%2C+6%7D%2C+%7B7%2C+8%7D%7D%7B%7B1%2C+2%7D%2C+%7B3%2C+4%7D%7D
    //  Verify the result of matrix multiplication
    Matrix2x2<int> ans2({{85, 126},
                         {193, 286}});
    REQUIRE(result2 == ans2);
}

TEST_CASE("Matrix2x2-Vector Multiplication") {
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});

    // Define a std::array vector of size 2
    std::array<int, 2> vec = {5, 6};

    // Perform matrix-vector multiplication
    std::array<int, 2> result = mat * vec;

    // Define the expected result
    //https://www.wolframalpha.com/input?i2d=true&i=%7B%7B1%2C2%7D%2C%7B3%2C4%7D%7D%7B%7B5%7D%2C%7B6%7D%7D
    std::array<int, 2> expected = {17, 39};  // Calculated as mat * vec

    // Check if the result matches the expected values
    REQUIRE(result == expected);
}

TEST_CASE("Matrix2x2    erminant")
{
    Matrix2x2<int> mat({{1, 2}, {3, 4}});

    // Calculate the determinant
    int det = mat.det();

    // Verify the determinant value
    REQUIRE(det == -2);
}

TEST_CASE("Matrix2x2 Addition")
{
    Matrix2x2<int> mat1({{1, 2},
                         {3, 4}});
    Matrix2x2<int> mat2({{5, 6},
                         {7, 8}});
    Matrix2x2<int> result = mat1 + mat2;

    Matrix2x2<int> ans({{6, 8},
                        {10, 12}});
    // Verify the result of matrix addition
    REQUIRE(result == ans);
}

TEST_CASE("Matrix2x2 Transpose")
{
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> transposed = mat.transpose();

    Matrix2x2<int> ans({{1, 3},
                        {2, 4}});
    // Verify the transpose of the matrix
    REQUIRE(transposed == ans);
}

TEST_CASE("Matrix2x2 Inverse")
{
    Matrix2x2<double> mat({{2.0, 1.0},
                           {1.0, 3.0}});
    Matrix2x2<double> inverse = mat.inverse();

    // Verify the inverse of the matrix
    REQUIRE(inverse[0][0] == doctest::Approx(3 / 5.0));
    REQUIRE(inverse[0][1] == doctest::Approx(-1 / 5.0));
    REQUIRE(inverse[1][0] == doctest::Approx(-1 / 5.0));
    REQUIRE(inverse[1][1] == doctest::Approx(2 / 5.0));
}
// Lagrange reduction matrixies
TEST_CASE("Matrix2x2 lag_m1 Test")
{
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> ans({{1, -2},
                        {3, -4}});
    mat.lag_m1();

    REQUIRE(mat == ans);
}

TEST_CASE("Matrix2x2 lag_m2 Test")
{
    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> ans({{2, 1},
                        {4, 3}});
    mat.lag_m2();

    REQUIRE(mat == ans);
}

TEST_CASE("Matrix2x2 lag_m3 Test")
{
    // https://www.wolframalpha.com/input?i=%7B%7B1%2C+2%7D%2C+%7B3%2C+4%7D%7D+%7B%7B1%2C+-1%7D%2C+%7B0%2C+1%7D%7D

    Matrix2x2<int> mat({{1, 2},
                        {3, 4}});
    Matrix2x2<int> ans({{1, 1},
                        {3, 1}});
    mat.lag_m3();

    REQUIRE(mat == ans);
}

TEST_CASE("Test similarityTransform function")
{
    Matrix2x2<double> a = {{1, 2},
                           {3, 4}};

    Matrix2x2<double> g = {{0.5, 0.2},
                           {-0.3, 0.8}};

    Matrix2x2<double> result = a.similarityTransform(g);
    // https://www.wolframalpha.com/input?i=inverse+%7B%7B0.5%2C+0.2%7D%2C+%7B-0.3%2C+0.8%7D%7D++%7B%7B1%2C2%7D%2C%7B3%2C4%7D%7D++%7B%7B0.5%2C+0.2%7D%2C+%7B-0.3%2C+0.8%7D%7D+

    double expected00 = -0.304348;
    double expected01 = 1.47826;
    double expected10 = 0.26087;
    double expected11 = 5.30435;

    REQUIRE(result[0][0] == doctest::Approx(expected00));
    REQUIRE(result[0][1] == doctest::Approx(expected01));
    REQUIRE(result[1][0] == doctest::Approx(expected10));
    REQUIRE(result[1][1] == doctest::Approx(expected11));
}

TEST_CASE("Test congruenceTransform function")
{
    Matrix2x2<double> a = {{1, 2},
                           {3, 4}};

    Matrix2x2<double> g = Matrix2x2<double>::rotation_matrix(3.0 / 4.0 * M_PI);

    Matrix2x2<double> result = a.congruenceTransform(g);
    // https://www.wolframalpha.com/input?i=transpose+%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D++%7B%7B1%2C2%7D%2C%7B3%2C4%7D%7D+++%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D+where+a%3D3pi%2F4+
    double expected00 = 0;
    double expected01 = -2;
    double expected10 = -1;
    double expected11 = 5;

    REQUIRE(result[0][0] == doctest::Approx(expected00));
    REQUIRE(result[0][1] == doctest::Approx(expected01));
    REQUIRE(result[1][0] == doctest::Approx(expected10));
    REQUIRE(result[1][1] == doctest::Approx(expected11));
}

TEST_CASE("Test isOrthogonal function")
{
    // Create an orthogonal matrix and test 'isOrthogonal' function.
    Matrix2x2<double> g = {{1, 0},
                           {0, -1}};

    REQUIRE(g.isOrthogonal() == true);

    g = Matrix2x2<double>(3);
    REQUIRE(g.isOrthogonal() == false);

    // 3 is a random theta. No special meaning.
    g = Matrix2x2<double>::reflection_matrix(2.0);
    REQUIRE(g.isOrthogonal() == true);
}

TEST_CASE("Test sym_orth_conjugate function")
{
    // Create a symmetric matrix and an orthogonal matrix and test 'sym_orth_conjugate' function.
    // Make sure to define your expected result based on the known computation of the symmetric orthogonal similarityTransform.
    // Replace the following values with the expected values.
    Matrix2x2<double> a = {{1, 2},
                           {2, 3}};

    Matrix2x2<double> g = Matrix2x2<double>::rotation_matrix(3.0 / 4.0 * M_PI);
    Matrix2x2<double> result = a.sym_orth_conjugate(g);

    // https://www.wolframalpha.com/input?i=inverse+%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D++%7B%7B1%2C2%7D%2C%7B2%2C3%7D%7D+++%7B%7Bcos%28a%29%2C+-sin%28a%29%7D%2C+%7Bsin%28a%29%2C+cos%28a%29%7D%7D+where+a%3D3pi%2F4+
    double expected00 = 0;
    double expected01 = -1;
    double expected10 = -1;
    double expected11 = 4;

    REQUIRE(result[0][0] == doctest::Approx(expected00));
    REQUIRE(result[0][1] == doctest::Approx(expected01));
    REQUIRE(result[1][0] == doctest::Approx(expected10));
    REQUIRE(result[1][1] == doctest::Approx(expected11));
}
