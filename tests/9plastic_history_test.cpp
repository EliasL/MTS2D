#include "Data/cereal_help.h"
#include "Mesh/mesh.h"
#include "Simulation/simulation.h"
#include "element_history_export.h"
#include "run/doctest.h"
#include <filesystem>
#include <fstream>
#include <sstream>

TEST_CASE("Edge flips inherit stored plastic history under either distance choice") {
  for (bool minimum : {false, true}) {
    CAPTURE(minimum);
    Mesh mesh(2, 2, false, "minor");
    mesh.minimum_transport_distance = minimum;
    mesh.applyTransformation(getShear(1));
    mesh.ensureFull();

    Matrix2d vertical = getShear(2).transpose();
    mesh.elements[0].T_p = vertical * getShear(1);
    mesh.elements[1].T_p = getShear(-1) * vertical;
    const Mesh before = mesh;
    mesh.flipEdge(mesh.elements[0], mesh.elements[1]);
    mesh.markDirty();
    mesh.ensureFull();

    CHECK(mesh.rmsDistanceToMesh(before) == 0);
    bool plasticFactorChanged = false;
    for (size_t i = 0; i < mesh.elements.size(); ++i) {
      CHECK(mesh.elements[i].T_p == before.elements[i].T_p);
      CHECK(mesh.elements[i].totalDeformation().isApprox(
          mesh.elements[i].F_E * before.elements[i].T_p, 1e-12));
      plasticFactorChanged |= mesh.elements[i].F_P != before.elements[i].F_P;
    }
    CHECK(plasticFactorChanged);

    const Matrix2d history = mesh.elements[0].T_p;
    mesh.ensureFull();
    mesh.markDirty();
    mesh.ensureFull();
    CHECK(mesh.elements[0].T_p == history);
  }
}

TEST_CASE("Force-only dumps preserve the pending plastic history update") {
  Mesh mesh(2, 2, false, "minor");
  mesh.applyTransformation(getShear(1));
  mesh.ensureGeometry();
  mesh.flipEdge(mesh.elements[0], mesh.elements[1]);
  mesh.markDirty();
  mesh.ensureFull();

  const Matrix2d previousHistory = mesh.elements[0].T_p;
  const Matrix2d previousPlastic = mesh.elements[0].F_P;
  mesh.applyTransformation(getShear(2).transpose());
  mesh.ensureForces();
  REQUIRE(mesh.getUpdateState() == Mesh::UpdateState::Forces);
  const std::string xml = serializeToXml(mesh);
  CHECK(mesh.getUpdateState() == Mesh::UpdateState::Forces);
  CHECK(mesh.elements[0].T_p == previousHistory);
  CHECK(mesh.elements[0].F_P == previousPlastic);

  Mesh loaded;
  deserializeFromXml(xml, loaded);
  mesh.ensureFull();
  loaded.ensureFull();
  CHECK(mesh.elements[0].F_P != previousPlastic);
  for (size_t i = 0; i < mesh.elements.size(); ++i) {
    CHECK(loaded.elements[i].T_p.isApprox(mesh.elements[i].T_p, 1e-12));
    CHECK(loaded.elements[i].totalDeformation().isApprox(
        mesh.elements[i].totalDeformation(), 1e-12));
  }
  CHECK(xml.find("<T_p>") != std::string::npos);
  CHECK(xml.find("F_P_H") == std::string::npos);
}

TEST_CASE("Setting a new reference resets accumulated plastic history") {
  Mesh mesh(3, 3, true, "major");
  mesh.applyTransformation(getShear(2));
  mesh.ensureFull();
  REQUIRE(mesh.elements[0].T_p != Matrix2d::Identity());
  mesh.setRefConfiguration();
  mesh.ensureFull();
  for (const auto &element : mesh.elements) {
    CHECK(element.T_p.isIdentity(1e-12));
    CHECK(element.totalDeformation().isIdentity(1e-12));
  }
}

TEST_CASE("Matrix history reports use the stored T_p and history-corrected T") {
  Mesh mesh(2, 2, false, "minor");
  mesh.applyTransformation(getShear(1));
  mesh.ensureFull();
  ElementStepReconnectSnapshot snapshot;
  snapshot.elementChanged = true;
  snapshot.before = captureElementMatrixHistoryRow(mesh, 0);
  mesh.flipEdge(mesh.elements[0], mesh.elements[1]);
  mesh.markDirty();
  mesh.ensureFull();
  snapshot.after = captureElementMatrixHistoryRow(mesh, 0);

  Config config;
  config.setDefaultValues();
  config.rows = config.cols = 2;
  const auto path = std::filesystem::path("test_data/plasticHistory/report.json");
  std::filesystem::create_directories(path.parent_path());
  writeElementStepReconnectJson(path, 0, "Plastic history across an edge flip",
                                {snapshot}, config);
  std::ifstream input(path);
  std::ostringstream contents;
  contents << input.rdbuf();
  CHECK(contents.str().find("\"T_p\":") != std::string::npos);
  CHECK(contents.str().find("\"T\":") != std::string::npos);
  CHECK(contents.str().find("\"H\":") == std::string::npos);
}
