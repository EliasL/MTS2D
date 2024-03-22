#include "run/doctest.h"
#include "../src/Mesh/node.h" // Include the header for your surface struct

// Test case for the NodeId struct
TEST_CASE("NodeId Struct Test")
{
    /*
    Visualization of col, row and i in NodeId for 3x3 matrix

    2       6   7   8
    1       3   4   5
    0       0   1   2

    y/x     0   1   2
    */
    NodeId id1 = {1, 1, 3};
    CHECK(id1.i == 4);
    NodeId id2 = {4, 3};
    CHECK(id2.col == 1);
    CHECK(id2.row == 1);
}

TEST_CASE("Node test")
{
    /*
    This is a node at the corner of a 2x2 mesh
    2 3
    0 1

    so it should have index i=3
    */
    Node n(1, 1, 1, 2);
    CHECK(n.id.i == 3);
    CHECK(n.id.col == 1);
    CHECK(n.id.row == 1);
    CHECK(n.pos()[0] == 1);
    CHECK(n.pos()[1] == 1);
    CHECK(n.init_pos()[0] == 1);
    CHECK(n.init_pos()[1] == 1);
}