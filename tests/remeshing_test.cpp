#include "../src/Data/data_export.h"
#include "../src/Mesh/mesh.h" // Include the header for your surface struct
#include "Mesh/compare_macros.h"
#include "Mesh/node.h"
#include "run/doctest.h"
#include <iostream>
#include <string>

void save(Mesh &mesh, std::string name) {
  mesh.updateMesh();
  writeMeshToVtu(mesh, "remeshing", "", name);
}

TEST_CASE("Mesh Initialization") {
  // Create a mesh
  Mesh mesh(2, 2, false);

  save(mesh, "InitialState");

  // Move upper nodes
  for (int i = 0; i < mesh.nodes.rows(); ++i) {
    for (int j = 0; j < mesh.nodes.cols(); ++j) {
      Node &n = mesh.nodes(i, j); // Access node reference
      if (i == 1) {
        n.setDisplacement({-0.3, 0});
      }
    }
  }
  save(mesh, "displaced");

  Mesh oldMesh = mesh;

  // Remesh
  mesh.setDiagonal(0, 0, false);

  save(mesh, "remeshed");

  // Compare the energy and forces on the nodes
  for (int i = 0; i < mesh.nrNodes; i++) {
    std::cout << mesh.nodes(i).f << '\n';
    CHECK(mesh.nodes(i).f == oldMesh.nodes(i).f);
  }
}
