#include "data_export.h"

#include "Data/lean_vtk.h"
#include "Mesh/mesh.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "settings.h"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// Get file path, and check that the path exists
std::string getFilePath(std::string fileName, std::string folderName,
                        std::string dataPath, std::string subDataFolder = "",
                        std::string fileType = ".vtu") {
  std::string fileNameWithType = fileName + fileType;
  std::string directory = getDataPath(folderName, dataPath);
  if (!subDataFolder.empty()) {
    directory += "/" + subDataFolder;
    fs::create_directories(directory);
  }

  // Check if the directory exists
  if (!fs::exists(directory) || !fs::is_directory(directory)) {
    throw std::runtime_error("Directory does not exist: " + directory +
                             "\nHave you run the createDataFolder function?");
  }

  return fs::path(directory) / fileNameWithType;
}

// If we want to store some data that does not depend on either the node or
// cell, it is inefficient to store the data multiple times. The simplest way I
// have found to store extra data is by including it in the file name, dataPath.
// Example: The variable foo and bar are stored as "_foo=0.32_bar=4_".
std::string makeFileName(const Mesh &mesh, std::string name) {
  std::stringstream ss;
  ss << name << "_load=" << mesh.load
     << "_nrM=" << mesh.nr_elements_with_m3_change << '_';
  return ss.str();
}

struct Vector2iHash {
  std::size_t operator()(const Vector2i &v) const noexcept {
    return std::hash<int>()(v.x()) ^ (std::hash<int>()(v.y()) << 1);
  }
};

// Gives each unique ghostNode an index
std::unordered_map<Vector2i, int, Vector2iHash>
constructGhostNodeIndexes(const Mesh &mesh) {
  // This is the maximum number of nodes that COULD be used to represent
  // a periodic mesh (within reason, technically, it's 3*nrElements)
  int maxNrNodes = mesh.cols * mesh.rows + 2 * (mesh.cols + mesh.rows);
  int currentIndex = 0;

  // Here we make sure we only add unique ghost nodes (uniqueness is checked by
  // row, col and periodic shift)
  std::unordered_map<Vector2i, int, Vector2iHash> uniqueGhostNodes(maxNrNodes);
  for (const TElement &e : mesh.elements) {
    for (const GhostNode &gn : e.ghostNodes) {
      // If the map does not yet contain this ghost node, we add it to the list
      // and give it an index
      if (uniqueGhostNodes.find(gn.id) == uniqueGhostNodes.end()) {
        uniqueGhostNodes[gn.id] = currentIndex;
        currentIndex++;
      }
    }
  }
  return uniqueGhostNodes;
}

namespace {
const char *vtuFieldLevelSuffix(VtuFieldLevel level) {
  switch (level) {
  case VtuFieldLevel::Minimal:
    return "minimal";
  case VtuFieldLevel::Extras:
    return "extras";
  case VtuFieldLevel::All:
  default:
    return "";
  }
}

bool sameExportPosition(const Vector2d &a, const Vector2d &b,
                        double tol = 1e-10) {
  return (a - b).cwiseAbs().maxCoeff() <= tol;
}

void validateRepeatedGhostNodeForExport(
    const GhostNode &gn, int nodeIndex,
    const std::vector<Vector2d> &storedPositions,
    const std::vector<int> &storedReferenceIds,
    const std::vector<Vector2i> &storedPeriodicShifts) {
  if (storedReferenceIds[nodeIndex] != gn.referenceId.i) {
    throw std::runtime_error(
        "writeMeshToVtu: repeated ghost node id maps to different "
        "reference nodes.");
  }
  if (storedPeriodicShifts[nodeIndex] != gn.periodicShift) {
    throw std::runtime_error(
        "writeMeshToVtu: repeated ghost node id maps to different periodic "
        "shifts.");
  }
  if (!sameExportPosition(storedPositions[nodeIndex], gn.pos)) {
    throw std::runtime_error(
        "writeMeshToVtu: repeated ghost node id maps to different positions.");
  }
}
} // namespace

std::string writeMeshToVtu(const Mesh &mesh, std::string folderName,
                           std::string dataPath, std::string fileName,
                           bool minimizationStep, VtuFieldLevel level,
                           std::string nameSuffix, std::string subDataFolder,
                           bool useReferenceElements) {

  const int dim = 3;
  const int cell_size = 3;
  int nrElements = mesh.nrElements;
  if (mesh.F_P_H.size() != static_cast<size_t>(nrElements)) {
    throw std::runtime_error(
        "writeMeshToVtu: branch history size does not match number of "
        "elements.");
  }

  // Due to periodic reconnecting, the number of nodes we use to represent the
  // periodic elements is not constant. We therefore create a node list here
  std::unordered_map<Vector2i, int, Vector2iHash> nodeMap;
  int nrNodes = 0;
  if (useReferenceElements) {
    nrNodes = nrElements * cell_size;
  } else {
    nodeMap = constructGhostNodeIndexes(mesh);
    nrNodes = nodeMap.size();
  }
  std::string subFolder = subDataFolder;
  if (fileName == "") {
    fileName = makeFileName(mesh, folderName);
  }

  if (level != VtuFieldLevel::All) {
    if (!fileName.empty() && fileName.back() != '_') {
      fileName += "_";
    }
    fileName += std::string(vtuFieldLevelSuffix(level)) + "_";
  }
  if (minimizationStep) {
    // Here we an extra folder for each minimization step
    if (!subFolder.empty()) {
      subFolder += "/";
    }
    subFolder += getMinDataSubFolder(mesh);
    fileName += "minStep=" + std::to_string(mesh.nrMinItterations) + '.' +
                std::to_string(mesh.nrMinFunctionCalls);
  }

  if (!nameSuffix.empty()) {
    if (!fileName.empty() && fileName.back() != '_') {
      fileName += "_";
    }
    fileName += nameSuffix;
  }

  std::string filePath;

  filePath = getFilePath(fileName + "." + std::to_string(mesh.loadSteps),
                         folderName, dataPath, subFolder);
  // std::cout << filePath << '\n';

  std::vector<double> points(nrNodes * dim);
  const bool minimal = (level == VtuFieldLevel::Minimal);
  const bool all = (level == VtuFieldLevel::All);

  const bool wantMatrices = all;
  const bool wantSigmaMatrix = all;
  const bool wantQuadrants = all;
  const bool wantSigma12Only = minimal;
  const bool wantFixed = !minimal;
  const bool wantRefIndex = !minimal;
  const bool wantDet = !minimal;
  const bool wantAngles = !minimal;
  const bool wantDisplacement = all;

  std::vector<double> force(nrNodes * dim);
  std::vector<double> displacement;
  // boolean values represented by 0.0 and 1.0
  std::vector<double> fixed;
  std::vector<double> refIndex; // Index of reference node
  if (wantFixed) {
    fixed.resize(nrNodes);
  }
  if (wantRefIndex) {
    refIndex.resize(nrNodes);
  }
  if (wantDisplacement) {
    displacement.resize(nrNodes * dim);
  }
  std::vector<int> connectivity(nrElements * cell_size); // Mesh order
  std::vector<double> energy(nrElements);
  std::vector<double> det;
  std::vector<double> F11, F12, F21, F22;
  std::vector<double> F_P11, F_P12, F_P21, F_P22;
  std::vector<double> F_E11, F_E12, F_E21, F_E22;
  std::vector<double> C11, C12, C22;
  std::vector<double> G11, G12, G22;
  std::vector<double> P11, P12, P21, P22;
  std::vector<double> T11, T12, T21, T22;
  std::vector<double> sigma11, sigma12, sigma22;

  std::vector<double> largeAngle;
  std::vector<double> smallAngle;

  // These should be int, but the library i am using only takes doubles
  std::vector<double> nrm3;               // Int
  std::vector<double> deltaNrm3;          // Int
  std::vector<double> twinID;             // Int
  std::vector<double> m11, m12, m21, m22; // Int
  std::vector<double> redQuadrant;        // Int

  if (wantDet) {
    det.resize(nrElements);
  }
  if (wantMatrices) {
    F11.resize(nrElements);
    F12.resize(nrElements);
    F21.resize(nrElements);
    F22.resize(nrElements);
    F_P11.resize(nrElements);
    F_P12.resize(nrElements);
    F_P21.resize(nrElements);
    F_P22.resize(nrElements);
    F_E11.resize(nrElements);
    F_E12.resize(nrElements);
    F_E21.resize(nrElements);
    F_E22.resize(nrElements);
    C11.resize(nrElements);
    C12.resize(nrElements);
    C22.resize(nrElements);
    G11.resize(nrElements);
    G12.resize(nrElements);
    G22.resize(nrElements);
    P11.resize(nrElements);
    P12.resize(nrElements);
    P21.resize(nrElements);
    P22.resize(nrElements);
    T11.resize(nrElements);
    T12.resize(nrElements);
    T21.resize(nrElements);
    T22.resize(nrElements);
    m11.resize(nrElements);
    m12.resize(nrElements);
    m21.resize(nrElements);
    m22.resize(nrElements);
  }
  if (wantSigmaMatrix) {
    sigma11.resize(nrElements);
    sigma12.resize(nrElements);
    sigma22.resize(nrElements);
  } else if (wantSigma12Only) {
    sigma12.resize(nrElements);
  }
  if (wantAngles) {
    largeAngle.resize(nrElements);
    smallAngle.resize(nrElements);
  }
  nrm3.resize(nrElements);
  deltaNrm3.resize(nrElements);
  twinID.resize(nrElements);
  if (wantQuadrants) {
    redQuadrant.resize(nrElements);
  }

  leanvtk::VTUWriter writer;

  // Instead of getting the data directly from the nodes in the mesh, we extract
  // the data from the nodes in the elements in the mesh. This is because they
  // have a displaced position and to not result in overlapping elements
  // DO NOT USE std::vector<bool>! This leads to memory
  // coruption errors that are difficult to track down
  std::vector<char> alreadyCopied(nrNodes, false);
  std::vector<Vector2d> storedPositions(nrNodes, Vector2d::Zero());
  std::vector<int> storedReferenceIds(nrNodes, -1);
  std::vector<Vector2i> storedPeriodicShifts(nrNodes, Vector2i::Zero());
  // Iterate over each element in the mesh
  for (int elementIndex = 0; elementIndex < nrElements; ++elementIndex) {
    const TElement &e = mesh.elements[elementIndex];
    const Vector2d referenceCentroidShift =
        useReferenceElements ? e.referenceCentroidShiftToCurrent()
                             : Vector2d::Zero();
    // Iterate over each node in the element
    for (size_t j = 0; j < e.ghostNodes.size(); ++j) {
      const GhostNode &gn = e.ghostNodes[j];
      int nodeIndex = 0;
      if (useReferenceElements) {
        nodeIndex = elementIndex * cell_size + static_cast<int>(j);
        // For the reference VTU, we only change the visualization gauge:
        // translate each disconnected reference triangle so its centroid
        // matches the current triangle centroid. This does not mutate the
        // stored ghost-node reference positions used by the simulation.
        const Vector2d exportedRefPos = gn.ref_pos + referenceCentroidShift;
        points[nodeIndex * dim + 0] = exportedRefPos[0];
        points[nodeIndex * dim + 1] = exportedRefPos[1];
        points[nodeIndex * dim + 2] = 0;
        force[nodeIndex * dim + 0] = gn.f[0];
        force[nodeIndex * dim + 1] = gn.f[1];
        force[nodeIndex * dim + 2] = 0;
        if (wantDisplacement) {
          const Vector2d exportedU = gn.pos - exportedRefPos;
          displacement[nodeIndex * dim + 0] = exportedU[0];
          displacement[nodeIndex * dim + 1] = exportedU[1];
          displacement[nodeIndex * dim + 2] = 0;
        }
        if (wantFixed) {
          fixed[nodeIndex] = mesh[gn]->fixedNode;
        }
        if (wantRefIndex) {
          refIndex[nodeIndex] = gn.referenceId.i;
        }
      } else {
        // Element index
        auto it = nodeMap.find(gn.id);
        assert(it != nodeMap.end()); // or handle the error
        if (it == nodeMap.end()) {
          throw std::runtime_error("Ghost node not found in node map.");
        }
        nodeIndex = it->second;
        if (!alreadyCopied[nodeIndex]) {
          storedPositions[nodeIndex] = gn.pos;
          storedReferenceIds[nodeIndex] = gn.referenceId.i;
          storedPeriodicShifts[nodeIndex] = gn.periodicShift;
          points[nodeIndex * dim + 0] = gn.pos[0];
          points[nodeIndex * dim + 1] = gn.pos[1];
          points[nodeIndex * dim + 2] = 0;
          // We need to take the force from the mesh nodes, not the element nodes
          const Node &n = *mesh[gn];
          force[nodeIndex * dim + 0] = n.f[0];

          force[nodeIndex * dim + 1] = n.f[1];
          force[nodeIndex * dim + 2] = 0;
          if (wantDisplacement) {
            const Vector2d &u = n.u();
            displacement[nodeIndex * dim + 0] = u[0];
            displacement[nodeIndex * dim + 1] = u[1];
            displacement[nodeIndex * dim + 2] = 0;
          }
          if (wantFixed) {
            fixed[nodeIndex] = n.fixedNode;
          }
          if (wantRefIndex) {
            refIndex[nodeIndex] = gn.referenceId.i;
          }
          alreadyCopied[nodeIndex] = true;
        } else {
          validateRepeatedGhostNodeForExport(gn, nodeIndex, storedPositions,
                                             storedReferenceIds,
                                             storedPeriodicShifts);
        }
      }
      // Here we define the meshing
      connectivity[elementIndex * cell_size + j] = nodeIndex;
    }

    energy[elementIndex] = e.energy;
    if (wantDet) {
      det[elementIndex] = e.C.determinant();
    }

    if (wantMatrices) {
      const Matrix2d T =
          e.totalBranch(mesh.F_P_H[static_cast<size_t>(elementIndex)]);
      F11[elementIndex] = e.F(0, 0);
      F12[elementIndex] = e.F(0, 1);
      F21[elementIndex] = e.F(1, 0);
      F22[elementIndex] = e.F(1, 1);
      F_P11[elementIndex] = e.F_P(0, 0);
      F_P12[elementIndex] = e.F_P(0, 1);
      F_P21[elementIndex] = e.F_P(1, 0);
      F_P22[elementIndex] = e.F_P(1, 1);
      F_E11[elementIndex] = e.F_E(0, 0);
      F_E12[elementIndex] = e.F_E(0, 1);
      F_E21[elementIndex] = e.F_E(1, 0);
      F_E22[elementIndex] = e.F_E(1, 1);

      C11[elementIndex] = e.C(0, 0);
      C12[elementIndex] = e.C(0, 1);
      C22[elementIndex] = e.C(1, 1);

      G11[elementIndex] = e.G(0, 0);
      G12[elementIndex] = e.G(0, 1);
      G22[elementIndex] = e.G(1, 1);
      P11[elementIndex] = e.P(0, 0);
      P12[elementIndex] = e.P(0, 1);
      P21[elementIndex] = e.P(1, 0);
      P22[elementIndex] = e.P(1, 1);
      T11[elementIndex] = T(0, 0);
      T12[elementIndex] = T(0, 1);
      T21[elementIndex] = T(1, 0);
      T22[elementIndex] = T(1, 1);
      m11[elementIndex] = e.M_l(0, 0);
      m12[elementIndex] = e.M_l(0, 1);
      m21[elementIndex] = e.M_l(1, 0);
      m22[elementIndex] = e.M_l(1, 1);
    }

    if (wantSigmaMatrix) {
      sigma11[elementIndex] = e.sigma(0, 0);
      sigma12[elementIndex] = e.sigma(0, 1);
      sigma22[elementIndex] = e.sigma(1, 1);
    } else if (wantSigma12Only) {
      sigma12[elementIndex] = e.sigma(0, 1);
    }
    if (wantAngles) {
      largeAngle[elementIndex] = e.largestAngle;
      smallAngle[elementIndex] = e.smallestAngle;
    }
    nrm3[elementIndex] = e.m3Nr;
    deltaNrm3[elementIndex] = e.m3Nr - e.pastM3Nr;
    twinID[elementIndex] = e.getElementTwin(mesh);
    if (wantQuadrants) {
      redQuadrant[elementIndex] = e.red_quadrant;
    }
  }

  // Debug confirm that all the nodes have been written to
  if (!useReferenceElements) {
    assert(std::all_of(alreadyCopied.begin(), alreadyCopied.end(),
                       [](bool value) { return value; }));
  }

  std::map<std::array<int, 3>, int> seenCells;
  if (!useReferenceElements) {
    for (int elementIndex = 0; elementIndex < nrElements; ++elementIndex) {
      std::array<int, 3> cell = {
          connectivity[elementIndex * cell_size + 0],
          connectivity[elementIndex * cell_size + 1],
          connectivity[elementIndex * cell_size + 2]};
      std::sort(cell.begin(), cell.end());
      auto [it, inserted] = seenCells.emplace(cell, elementIndex);
      if (!inserted) {
        throw std::runtime_error(
            "writeMeshToVtu: duplicate exported cell connectivity for "
            "elements " +
            std::to_string(it->second) + " and " +
            std::to_string(elementIndex) + ".");
      }
    }
  }

  // connect data to writer
  writer.add_cell_scalar_field("energy_field", energy);
  if (wantDet) {
    writer.add_cell_scalar_field("det", det);
  }
  if (wantFixed) {
    writer.add_scalar_field("fixed", fixed);
  }
  if (wantRefIndex) {
    writer.add_scalar_field("refIndex", refIndex);
  }
  if (wantMatrices) {
    writer.add_cell_scalar_field("F11", F11);
    writer.add_cell_scalar_field("F12", F12);
    writer.add_cell_scalar_field("F21", F21);
    writer.add_cell_scalar_field("F22", F22);
    writer.add_cell_scalar_field("F_P11", F_P11);
    writer.add_cell_scalar_field("F_P12", F_P12);
    writer.add_cell_scalar_field("F_P21", F_P21);
    writer.add_cell_scalar_field("F_P22", F_P22);
    writer.add_cell_scalar_field("F_E11", F_E11);
    writer.add_cell_scalar_field("F_E12", F_E12);
    writer.add_cell_scalar_field("F_E21", F_E21);
    writer.add_cell_scalar_field("F_E22", F_E22);
    writer.add_cell_scalar_field("C11", C11);
    writer.add_cell_scalar_field("C12", C12);
    writer.add_cell_scalar_field("C22", C22);
    writer.add_cell_scalar_field("G11", G11);
    writer.add_cell_scalar_field("G12", G12);
    writer.add_cell_scalar_field("G22", G22);
    writer.add_cell_scalar_field("P11", P11);
    writer.add_cell_scalar_field("P12", P12);
    writer.add_cell_scalar_field("P21", P21);
    writer.add_cell_scalar_field("P22", P22);
    writer.add_cell_scalar_field("T11", T11);
    writer.add_cell_scalar_field("T12", T12);
    writer.add_cell_scalar_field("T21", T21);
    writer.add_cell_scalar_field("T22", T22);
    writer.add_cell_scalar_field("m11", m11);
    writer.add_cell_scalar_field("m12", m12);
    writer.add_cell_scalar_field("m21", m21);
    writer.add_cell_scalar_field("m22", m22);
  }
  if (wantSigmaMatrix) {
    writer.add_cell_scalar_field("sigma11", sigma11);
    writer.add_cell_scalar_field("sigma12", sigma12);
    writer.add_cell_scalar_field("sigma22", sigma22);
  } else if (wantSigma12Only) {
    writer.add_cell_scalar_field("sigma12", sigma12);
  }
  if (wantAngles) {
    writer.add_cell_scalar_field("smallAngle", smallAngle);
    writer.add_cell_scalar_field("largeAngle", largeAngle);
  }

  writer.add_cell_scalar_field("nrm3", nrm3);
  writer.add_cell_scalar_field("deltaNrm3", deltaNrm3);
  writer.add_cell_scalar_field("twinID", twinID);
  if (wantQuadrants) {
    writer.add_cell_scalar_field("red_quadrant", redQuadrant);
  }

  writer.add_vector_field("stress_field", force, dim);
  if (wantDisplacement) {
    writer.add_vector_field("displacement", displacement, dim);
  }

  // write data
  writer.write_surface_mesh(filePath, dim, cell_size, points, connectivity);
  // Save a compressed version
  // Paraview doesn't open compressed files >:(
  // TODO, if used, remember to update the dataPath properly so pvd collection
  // is updated.
  // filePath = compressFile(filePath);
  return filePath;
}

void createCollection(const std::string &folderPath,
                      const std::string &destination,
                      const std::string &regexPattern,
                      const std::string &extension,
                      const std::vector<double> &timestep,
                      const std::string &collectionFileStem,
                      const std::string &requiredSubstring,
                      const std::string &excludedSubstring) {

  // Convert input paths to absolute to avoid ambiguities.
  fs::path absFolder = fs::absolute(folderPath);
  fs::path absDestination = fs::absolute(destination);

  // Check that the input folder exists.
  if (!fs::exists(absFolder)) {
    throw std::runtime_error("The folder path does not exist: " +
                             absFolder.string());
  }

  struct TimedFile {
    double baseTime = 0.0;
    fs::path path;
    fs::file_time_type writeTime;
  };

  // Container for storing files along with the extracted number from their
  // filenames and their write time.
  std::vector<TimedFile> timedFiles;
  // last number Example filename:
  // simpleShear,s50x50l0.297,1e-05,0.3PBCt6epsR1e-05LBFGSEpsg1e-08logDuringMinimization1energyDropThreshold1e-10s0_load=0.29701_nrM=0_minimal_minStep=432.7_post.14702
  std::regex regex;
  if (regexPattern.empty()) {
    regex = std::regex(".*\\.([0-9]+)\\.vtu");
  } else {
    regex = std::regex(regexPattern);
  }

  // Iterate over the directory entries.
  for (const auto &entry : fs::directory_iterator(absFolder)) {
    // Only process regular files with the given extension.
    if (entry.is_regular_file() && entry.path().extension() == extension) {
      const std::string filename = entry.path().filename().string();
      if (!requiredSubstring.empty() &&
          filename.find(requiredSubstring) == std::string::npos) {
        continue;
      }
      if (!excludedSubstring.empty() &&
          filename.find(excludedSubstring) != std::string::npos) {
        continue;
      }
      std::smatch match;
      if (std::regex_match(filename, match, regex) && match.size() == 2) {
        double number = std::stod(match[1].str());
        timedFiles.push_back(
            TimedFile{number, entry.path(), fs::last_write_time(entry.path())});
      }
    }
  }

  // Sort the files based on the extracted numeric value.
  std::sort(timedFiles.begin(), timedFiles.end(), [](const TimedFile &a,
                                                     const TimedFile &b) {
    if (a.baseTime != b.baseTime) {
      return a.baseTime < b.baseTime;
    }
    if (a.writeTime != b.writeTime) {
      return a.writeTime < b.writeTime;
    }
    return a.path.string() < b.path.string();
  });

  // Construct the full path of the collection file (output file).
  fs::path collectionFilePath = absDestination / (collectionFileStem + ".pvd");
  std::ofstream outFile(collectionFilePath);
  if (!outFile) {
    throw std::runtime_error("Failed to create output file: " +
                             collectionFilePath.string());
  }

  // Write the XML header and open the Collection tag.
  outFile << "<?xml version=\"1.0\"?>\n";
  outFile << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
  outFile << "  <Collection>\n";

  // Process each file and compute its relative path with respect to the
  // destination.
  std::size_t i = 0;
  while (i < timedFiles.size()) {
    const double baseTime =
        (i < timestep.size()) ? timestep[i] : timedFiles[i].baseTime;
    std::size_t groupEnd = i + 1;
    while (groupEnd < timedFiles.size()) {
      const double nextBaseTime =
          (groupEnd < timestep.size()) ? timestep[groupEnd]
                                       : timedFiles[groupEnd].baseTime;
      if (nextBaseTime != baseTime) {
        break;
      }
      ++groupEnd;
    }

    const std::size_t groupSize = groupEnd - i;
    const double duplicateWindow = 0.1;

    for (std::size_t groupIndex = 0; groupIndex < groupSize; ++groupIndex) {
      double timeValue = baseTime;
      if (groupSize > 1) {
        const double fraction =
            static_cast<double>(groupIndex) / static_cast<double>(groupSize - 1);
        timeValue += duplicateWindow * fraction;
      }

      // Compute the absolute file path and then its relative path from the
      // destination directory.
      fs::path absFilePath = fs::absolute(timedFiles[i + groupIndex].path);
      fs::path fileRelativePath = fs::relative(absFilePath, absDestination);

      // Write the DataSet element, ensuring the file attribute is quoted.
      outFile << "    <DataSet timestep=\"" << timeValue
              << "\" group=\"\" part=\"0\" file=\""
              << fileRelativePath.string() << "\"/>\n";
    }

    i = groupEnd;
  }

  // Close the XML tags.
  outFile << "  </Collection>\n";
  outFile << "</VTKFile>\n";

  outFile.close();
}

