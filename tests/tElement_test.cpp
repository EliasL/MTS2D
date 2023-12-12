#include "run/doctest.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct

TEST_CASE("Mesh Initialization")
{
    // Create a surface with known dimensions for testing
    Mesh mesh(2, 2);

    // There should now be two elements.
    REQUIRE(mesh.nrElements == 2);
}

TEST_CASE("Update current state")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[1];

    Matrix2x2<double> expectedState = {{-1, 0},
                                       {1, 1}};

    REQUIRE(e.currentState == expectedState);
    REQUIRE(e.currentState * e.invRefState == Matrix2x2<double>::identity());
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
    REQUIRE(e.F == shear2*shear);
}

TEST_CASE("Update metric tensor")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 4},
                               {0, 1}};
    // https://www.wolframalpha.com/input?i=transpose%28%7B%7B1%2C4%7D%2C%7B0%2C1%7D%7D%29.%7B%7B1%2C4%7D%2C%7B0%2C1%7D%7D
    Matrix2x2<double> ans =   {{1, 4},
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

    mesh.applyTransformation(shear);
    e.update();

    REQUIRE(e.C_ == Matrix2x2<double>::identity());
}


TEST_CASE("Update energy and reduced stress")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};

    e.update();

    REQUIRE(e.r_s[0][0] == doctest::Approx(0));
    REQUIRE(e.r_s[0][1] == doctest::Approx(0));
    REQUIRE(e.r_s[1][0] == doctest::Approx(0));
    REQUIRE(e.r_s[1][1] == doctest::Approx(0));

    REQUIRE(e.energy == doctest::Approx(3.9116));


    mesh.applyTransformation(shear);
    e.update();
     //TODO
}


TEST_CASE("Update Piola stress")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};

    e.update();

    REQUIRE(e.P[0][0] == doctest::Approx(0));
    REQUIRE(e.P[0][1] == doctest::Approx(0));
    REQUIRE(e.P[1][0] == doctest::Approx(0));
    REQUIRE(e.P[1][1] == doctest::Approx(0));

    mesh.applyTransformation(shear);
    e.update();
    //TODO
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
    //TODO
    
}

