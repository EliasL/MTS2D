#include "run/doctest.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct

TEST_CASE("Mesh Initialization")
{
    // Create a surface with known dimensions for testing
    Mesh mesh(2, 2);

    // There should now be two elements.
    REQUIRE(mesh.nrElements == 2);
}

TEST_CASE("Update deformation gradiant")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 4},
                               {0, 1}};

    mesh.applyTransformation(shear);
    e.update();
    REQUIRE(e.F == shear);

    Matrix2x2<double> shear2 = {{1, 0},
                                {3, 1}};

    mesh.applyTransformation(shear2);
    e.update();
    REQUIRE(e.F == shear2 * shear);
}

TEST_CASE("Update metric tensor")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 4},
                               {0, 1}};
    // https://www.wolframalpha.com/input?i=transpose%28%7B%7B1%2C4%7D%2C%7B0%2C1%7D%7D%29.%7B%7B1%2C4%7D%2C%7B0%2C1%7D%7D
    Matrix2x2<double> ans = {{1, 4},
                             {4, 17}};

    mesh.applyTransformation(shear);
    e.update();

    REQUIRE(e.C == ans);
}

TEST_CASE("Update reduced metric tensor")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 4},
                               {0, 1}};
    Matrix2x2<double> C_Ans = {{1, 0},
                               {0, 1}};
    Matrix2x2<double> mAns = {{1, -4},
                              {0, 1}};

    mesh.applyTransformation(shear);
    e.update();

    REQUIRE(e.C_ == C_Ans);
    REQUIRE(e.m == mAns);

    // Use the lagrange redution to test a more difficult example
    C_Ans = {{1.05565452, 0.52639976},
             {0.52639976, 1.20976767}};
    // https://www.wolframalpha.com/input?i=%7B%7B0%2C1%7D%2C%7B1%2C0%7D%7D.%7B%7B1%2C-1%7D%2C%7B0%2C1%7D%7D.%7B%7B1%2C-1%7D%2C%7B0%2C1%7D%7D.%7B%7B1%2C0%7D%2C%7B0%2C-1%7D%7D.%7B%7B0%2C1%7D%2C%7B1%2C0%7D%7D
    mAns = {{-1, 0},
            {2, 1}};
    e = TElement::lagrangeReduction(3.78912615, 1.20976767, 1.89313557);
    REQUIRE(approxEqual(e.C_, C_Ans));
    REQUIRE(e.m == mAns);
}

TEST_CASE("Update energy and reduced stress")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];

    e.update();

    REQUIRE(e.r_s[0][0] == doctest::Approx(0));
    REQUIRE(e.r_s[0][1] == doctest::Approx(0));
    REQUIRE(e.r_s[1][0] == doctest::Approx(0));
    REQUIRE(e.r_s[1][1] == doctest::Approx(0));

    REQUIRE(e.energy == doctest::Approx(3.9116));

    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};
    mesh.applyTransformation(shear);
    e.update();
    // TODO
}

TEST_CASE("Update Piola stress")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];

    e.update();

    REQUIRE(e.P[0][0] == doctest::Approx(0));
    REQUIRE(e.P[0][1] == doctest::Approx(0));
    REQUIRE(e.P[1][0] == doctest::Approx(0));
    REQUIRE(e.P[1][1] == doctest::Approx(0));

    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};
    mesh.applyTransformation(shear);
    e.update();
    // TODO
}

TEST_CASE("Apply forces on nodes")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};

    e.update();

    mesh.applyTransformation(shear);
    e.update();
    // TODO
}
