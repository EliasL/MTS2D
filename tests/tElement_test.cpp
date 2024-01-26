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
    // We use a mesh to initialize an element.
    Mesh mesh(2, 2);
    TElement e = mesh.elements[0];

    e.update();

    REQUIRE(e.r_s[0][0] == doctest::Approx(0));
    REQUIRE(e.r_s[0][1] == doctest::Approx(0));
    REQUIRE(e.r_s[1][0] == doctest::Approx(0));
    REQUIRE(e.r_s[1][1] == doctest::Approx(0));

    // Based on simulation from Umut sent by email on Jan 23, 2024, 10:04 AM
    REQUIRE(e.energy == doctest::Approx(3.91162));

    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};
    mesh.applyTransformation(shear);
    e.update();

    // Based on simulation from Umut sent by email on Jan 23, 2024, 10:22 AM
    REQUIRE(e.energy == doctest::Approx(4.00204));
}

TEST_CASE("Update Piola stress")
{
    // We use a mesh to initialize an element.
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
    // We use a mesh to initialize elements.
    Mesh mesh(3, 3);
    std::vector<TElement> &e = mesh.elements;
    Matrix2x2<double> shear = {{1, 0.5},
                               {0, 1}};



    // Element 0 Node Forces
    std::cout << "Ement 0:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 1 Node Forces
    std::cout << "Element 1:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 2 Node Forces
    std::cout << "Element 2:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 3 Node Forces
    std::cout << "Element 3:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 4 Node Forces
    std::cout << "Element 4:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 5 Node Forces
    std::cout << "Element 5:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 6 Node Forces
    std::cout << "Element 6:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    // Element 7 Node Forces
    std::cout << "Element 7:\n";
    std::cout << "P: "<< e[0].P<< "\n";

    for (size_t i = 0; i < e.size(); i++)
    {
        e[i].update();
        e[i].applyForcesOnNodes();
        std::cout << "element: "<<i<<" nodes: " << e[i].n1->id.i << ", "<< e[i].n2->id.i << ", "<< e[i].n3->id.i <<'\n';
    }

 
 // Element 0 Node Forces
    std::cout << "Element 0:\n";
    std::cout << "F11: "<< e[0].F[0][0] <<" F12: "<< e[0].F[0][1] <<" F21: "<< e[0].F[1][0] <<" F22: "<< e[0].F[1][1] <<"\n";
    // Element 1 Node Forces
    std::cout << "Element 1:\n";
    std::cout << "F11: "<< e[1].F[0][0] <<" F12: "<< e[1].F[0][1] <<" F21: "<< e[1].F[1][0] <<" F22: "<< e[1].F[1][1] <<"\n";
    // Element 2 Node Forces
    std::cout << "Element 2:\n";
    std::cout << "F11: "<< e[2].F[0][0] <<" F12: "<< e[2].F[0][1] <<" F21: "<< e[2].F[1][0] <<" F22: "<< e[2].F[1][1] <<"\n";
    // Element 3 Node Forces
    std::cout << "Element 3:\n";
    std::cout << "F11: "<< e[3].F[0][0] <<" F12: "<< e[3].F[0][1] <<" F21: "<< e[3].F[1][0] <<" F22: "<< e[3].F[1][1] <<"\n";
    // Element 4 Node Forces
    std::cout << "Element 4:\n";
    std::cout << "F11: "<< e[4].F[0][0] <<" F12: "<< e[4].F[0][1] <<" F21: "<< e[4].F[1][0] <<" F22: "<< e[4].F[1][1] <<"\n";
    // Element 5 Node Forces
    std::cout << "Element 5:\n";
    std::cout << "F11: "<< e[5].F[0][0] <<" F12: "<< e[5].F[0][1] <<" F21: "<< e[5].F[1][0] <<" F22: "<< e[5].F[1][1] <<"\n";
    // Element 6 Node Forces
    std::cout << "Element 6:\n";
    std::cout << "F11: "<< e[6].F[0][0] <<" F12: "<< e[6].F[0][1] <<" F21: "<< e[6].F[1][0] <<" F22: "<< e[6].F[1][1] <<"\n";
    // Element 7 Node Forces
    std::cout << "Element 7:\n";
    std::cout << "F11: "<< e[7].F[0][0] <<" F12: "<< e[7].F[0][1] <<" F21: "<< e[7].F[1][0] <<" F22: "<< e[7].F[1][1] <<"\n";
    
    
    // Element 0 Node Forces
    std::cout << "Ement 0:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 1 Node Forces
    std::cout << "Element 1:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 2 Node Forces
    std::cout << "Element 2:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 3 Node Forces
    std::cout << "Element 3:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 4 Node Forces
    std::cout << "Element 4:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 5 Node Forces
    std::cout << "Element 5:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 6 Node Forces
    std::cout << "Element 6:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";

    // Element 7 Node Forces
    std::cout << "Element 7:\n";
    std::cout << "invJ: "<< e[0].invJacobianRef<< "\n";





    mesh.applyTransformation(shear);
    for (size_t i = 0; i < e.size(); i++)
    {
        e[i].update();
        e[i].applyForcesOnNodes();
         std::cout << "element: "<<i<< " r_s: " << e[i].r_s <<'\n';
    }
       std::cout << "Element 1:\n";
    std::cout << "F11: "<< e[1].F[0][0] <<" F12: "<< e[1].F[0][1] <<" F21: "<< e[1].F[1][0] <<" F22: "<< e[1].F[1][1] <<"\n";






    // Element 0 Node Forces
    REQUIRE(e[0].n1->f_x == doctest::Approx(0.0462536));
    REQUIRE(e[0].n1->f_y == doctest::Approx(-0.0231268));
    REQUIRE(e[0].n2->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[0].n2->f_y == doctest::Approx(-0.0925071));
    REQUIRE(e[0].n3->f_x == doctest::Approx(0.0925071));
    REQUIRE(e[0].n3->f_y == doctest::Approx(0.0462536));

    // Element 1 Node Forces
    REQUIRE(e[1].n1->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[1].n1->f_y == doctest::Approx(-0.0925071));
    REQUIRE(e[1].n2->f_x == doctest::Approx(-0.0462536));
    REQUIRE(e[1].n2->f_y == doctest::Approx(-0.0693803));
    REQUIRE(e[1].n3->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[1].n3->f_y == doctest::Approx(0));

    // Element 2 Node Forces
    REQUIRE(e[2].n1->f_x == doctest::Approx(0.0925071));
    REQUIRE(e[2].n1->f_y == doctest::Approx(0.0462536));
    REQUIRE(e[2].n2->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[2].n2->f_y == doctest::Approx(0));
    REQUIRE(e[2].n3->f_x == doctest::Approx(0.0462536));
    REQUIRE(e[2].n3->f_y == doctest::Approx(0.0693803));

    // Element 3 Node Forces
    REQUIRE(e[3].n1->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[3].n1->f_y == doctest::Approx(0));
    REQUIRE(e[3].n2->f_x == doctest::Approx(-0.0925071));
    REQUIRE(e[3].n2->f_y == doctest::Approx(-0.0462536));
    REQUIRE(e[3].n3->f_x == doctest::Approx(0));
    REQUIRE(e[3].n3->f_y == doctest::Approx(0.0925071));

    // Element 4 Node Forces
    REQUIRE(e[4].n1->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[4].n1->f_y == doctest::Approx(0));
    REQUIRE(e[4].n2->f_x == doctest::Approx(0.0925071));
    REQUIRE(e[4].n2->f_y == doctest::Approx(0.0462536));
    REQUIRE(e[4].n3->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[4].n3->f_y == doctest::Approx(-0.0925071));

    // Element 5 Node Forces
    REQUIRE(e[5].n1->f_x == doctest::Approx(-0.0925071));
    REQUIRE(e[5].n1->f_y == doctest::Approx(-0.0462536));
    REQUIRE(e[5].n2->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[5].n2->f_y == doctest::Approx(0));
    REQUIRE(e[5].n3->f_x == doctest::Approx(-0.0462536));
    REQUIRE(e[5].n3->f_y == doctest::Approx(-0.0693803));

    // Element 6 Node Forces
    REQUIRE(e[6].n1->f_x == doctest::Approx(0));
    REQUIRE(e[6].n1->f_y == doctest::Approx(0.0925071));
    REQUIRE(e[6].n2->f_x == doctest::Approx(0.0462536));
    REQUIRE(e[6].n2->f_y == doctest::Approx(0.0693803));
    REQUIRE(e[6].n3->f_x == doctest::Approx(3.46945e-18).epsilon(0.01));
    REQUIRE(e[6].n3->f_y == doctest::Approx(0));

    // Element 7 Node Forces
    REQUIRE(e[7].n1->f_x == doctest::Approx(-0.0462536));
    REQUIRE(e[7].n1->f_y == doctest::Approx(0.0231268));
    REQUIRE(e[7].n2->f_x == doctest::Approx(0));
    REQUIRE(e[7].n2->f_y == doctest::Approx(0.0925071));
    REQUIRE(e[7].n3->f_x == doctest::Approx(-0.0925071));
    REQUIRE(e[7].n3->f_y == doctest::Approx(-0.0462536));
}