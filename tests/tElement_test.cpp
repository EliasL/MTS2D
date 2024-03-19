#include "run/doctest.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct

TEST_CASE("Mesh Initialization")
{
    // Create a surface with known dimensions for testing
    Mesh mesh(2, 2, false);

    // There should now be two elements.
    CHECK(mesh.nrElements == 2);
    // Create a periodic surface with known dimensions for testing
    mesh = Mesh(2, 2, true);

    // There should now be eight elements.
    CHECK(mesh.nrElements == 8);
}

TEST_CASE("Update deformation gradiant")
{
    // We use a mesh to initialize an element. (Not best practice for testing)
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];
    Matrix2x2<double> shear = {{1, 4},
                               {0, 1}};

    mesh.applyTransformation(shear);
    mesh.updateElements();
    CHECK(e.F == shear);

    Matrix2x2<double> shear2 = {{1, 0},
                                {3, 1}};

    mesh.applyTransformation(shear2);
    mesh.updateElements();
    CHECK(e.F == shear2 * shear);
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
    mesh.updateElements();

    CHECK(e.C == ans);
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
    mesh.updateElements();

    CHECK(e.C_ == C_Ans);
    CHECK(e.m == mAns);

    // Use the lagrange redution to test a more difficult example
    C_Ans = {{1.05565452, 0.52639976},
             {0.52639976, 1.20976767}};
    // https://www.wolframalpha.com/input?i=%7B%7B0%2C1%7D%2C%7B1%2C0%7D%7D.%7B%7B1%2C-1%7D%2C%7B0%2C1%7D%7D.%7B%7B1%2C-1%7D%2C%7B0%2C1%7D%7D.%7B%7B1%2C0%7D%2C%7B0%2C-1%7D%7D.%7B%7B0%2C1%7D%2C%7B1%2C0%7D%7D
    mAns = {{-1, 0},
            {2, 1}};
    // e = TElement::lagrangeReduction(3.78912615, 1.20976767, 1.89313557);
    // CHECK(approxEqual(e.C_, C_Ans));
    // CHECK(e.m == mAns);
}

TEST_CASE("Update energy and reduced stress")
{
    // We use a mesh to initialize an element.
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];

    mesh.updateElements();

    CHECK(e.r_s[0][0] == doctest::Approx(0));
    CHECK(e.r_s[0][1] == doctest::Approx(0));
    CHECK(e.r_s[1][0] == doctest::Approx(0));
    CHECK(e.r_s[1][1] == doctest::Approx(0));

    // Validated by Umut's code
    CHECK(e.energy == doctest::Approx(3.91162));

    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};
    mesh.applyTransformation(shear);
    mesh.updateElements();

    // Validated by Umut's code
    CHECK(e.energy == doctest::Approx(4.00204));
}

TEST_CASE("Update Piola stress")
{
    // We use a mesh to initialize an element.
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];

    mesh.updateElements();

    CHECK(e.P[0][0] == doctest::Approx(0));
    CHECK(e.P[0][1] == doctest::Approx(0));
    CHECK(e.P[1][0] == doctest::Approx(0));
    CHECK(e.P[1][1] == doctest::Approx(0));

    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};
    mesh.applyTransformation(shear);
    mesh.updateElements();

    // Validated by Umut's code
    CHECK(e.P[0][0] == doctest::Approx(-0.0462536));
    CHECK(e.P[0][1] == doctest::Approx(-3.47E-18));
    CHECK(e.P[1][0] == doctest::Approx(-0.0231268));
    CHECK(e.P[1][1] == doctest::Approx(0.0462536));
}

TEST_CASE("Apply forces on nodes at rest")
{
    // We use a mesh to initialize elements.
    Mesh mesh(3, 3, false);

    mesh.updateElements();
    mesh.applyForceFromElementsToNodes();

    for (size_t i = 0; i < mesh.nodes.data.size(); i++)
    {
        CHECK(mesh.nodes.data[i].f_x() == doctest::Approx(0));
        CHECK(mesh.nodes.data[i].f_y() == doctest::Approx(0));
    }
}

TEST_CASE("Apply forces on nodes")
{
    // We use a mesh to initialize elements.
    Mesh mesh(3, 3, false);
    auto &e = mesh.elements;
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};

    mesh.applyTransformation(shear);
    mesh.updateElements();
    mesh.applyForceFromElementsToNodes();

    // Validated by Umut's code
    CHECK(mesh.nodes.data[0].f_x() == doctest::Approx(0.0462536));
    CHECK(mesh.nodes.data[0].f_y() == doctest::Approx(-0.0231268));

    CHECK(mesh.nodes.data[1].f_x() == doctest::Approx(3.46945e-18));
    CHECK(mesh.nodes.data[1].f_y() == doctest::Approx(-0.0925071));

    CHECK(mesh.nodes.data[2].f_x() == doctest::Approx(-0.0462536));
    CHECK(mesh.nodes.data[2].f_y() == doctest::Approx(-0.0693803));

    CHECK(mesh.nodes.data[3].f_x() == doctest::Approx(0.0925071));
    CHECK(mesh.nodes.data[3].f_y() == doctest::Approx(0.0462536));

    CHECK(mesh.nodes.data[4].f_x() == doctest::Approx(3.46945e-18));
    CHECK(mesh.nodes.data[4].f_y() == doctest::Approx(0));

    CHECK(mesh.nodes.data[5].f_x() == doctest::Approx(-0.0925071));
    CHECK(mesh.nodes.data[5].f_y() == doctest::Approx(-0.0462536));

    CHECK(mesh.nodes.data[6].f_x() == doctest::Approx(0.0462536));
    CHECK(mesh.nodes.data[6].f_y() == doctest::Approx(0.0693803));

    CHECK(mesh.nodes.data[7].f_x() == doctest::Approx(0));
    CHECK(mesh.nodes.data[7].f_y() == doctest::Approx(0.0925071));

    CHECK(mesh.nodes.data[8].f_x() == doctest::Approx(-0.0462536));
    CHECK(mesh.nodes.data[8].f_y() == doctest::Approx(0.0231268));
}