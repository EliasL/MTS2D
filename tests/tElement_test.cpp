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