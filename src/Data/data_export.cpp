#include "data_export.h"

#include "../Simulation/simulation.h"
#include "Data/lean_vtk.h"
#include "Data/param_parser.h"
#include "Eigen/Core"
#include "Mesh/mesh.h"
#include "Mesh/node.h"
#include "Mesh/tElement.h"
#include "settings.h"
#include <Eigen/Dense>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem; // Alias for filesystem

std::string findOutputPath() {
  // Define the paths to check
  std::vector<std::string> paths = {
      "/Volumes/data/",
      "/media/elias/dataStorage/",
      "/data/elundheim/",  // PMMH
      "/data2/elundheim/", // PMMH
      "/Users/elias/Work/PhD/Code/localData/",
      "/lustre/fswork/projects/rech/bph/uog82gz/", // JeanZay
  };

  // Initialize a variable to store the chosen path
  fs::path outputFolder = OUTPUTFOLDERPATH;
  fs::path pre_chosen;
  fs::path chosen;

  for (const auto &base : paths) {
    if (fs::exists(base / outputFolder)) {
      chosen = base;
      break;
    }
    if (pre_chosen.empty() && fs::exists(base)) {
      pre_chosen = base;
    }
  }

  if (chosen.empty()) {
    chosen = pre_chosen;
  }

  if (chosen.empty()) {
    throw std::runtime_error(
        "Out path does not exist. Is your storage device connected?");
  }

  chosen /= outputFolder;
  return chosen.string();
}

std::string searchForConfig(std::string dumpPath) {
  // Extract the dataFolder path
  std::string dataFolderPath =
      dumpPath.substr(0, dumpPath.rfind(DUMPFOLDERPATH));

  // Append config.conf to the dataFolder path
  std::string configFilePath = dataFolderPath + CONFIGNAME;

  // Check if the config.conf file exists
  if (fs::exists(configFilePath)) {
    return configFilePath;
  } else {
    return "";
  }
}

std::string getCurrentDate() {
  auto now = std::chrono::system_clock::now();
  std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
  std::tm now_tm = *std::localtime(&now_time_t);

  char buffer[20]; // Adjust size as needed for the format
  // This format uses 16 chars, so 17 chars would in theory be the
  // minimum buffer size needed, but we have enough memory, so 20 is safer.
  strftime(buffer, sizeof(buffer), "%H.%M~%d.%m.%Y", &now_tm);

  std::stringstream ss;
  ss << buffer;
  return ss.str();
}

std::string getFolderPath(const std::string &name, const std::string &dataPath,
                          const fs::path &subfolder = "") {
  // Construct the full path
  fs::path fullPath = fs::path(dataPath) / name / subfolder;

  // Ensure that the directories exist
  fs::create_directories(fullPath);

  // Convert to string and return with trailing slash
  return fullPath.string();
}

std::string getOutputPath(const std::string &name,
                          const std::string &dataPath) {
  return getFolderPath(name, dataPath);
}

std::string getDataPath(const std::string &name, const std::string &dataPath) {
  return getFolderPath(name, dataPath, DATAFOLDERPATH);
}
std::string getMinDataSubFolder(const Mesh &mesh) {

  return std::string(MINDATAFOLDERPATH) + "/step" +
         std::to_string(mesh.loadSteps);
}

std::string getFramePath(const std::string &name, const std::string &dataPath) {
  return getFolderPath(name, dataPath, FRAMEFOLDERPATH);
}

std::string getDumpPath(const std::string &name, const std::string &dataPath) {
  return getFolderPath(name, dataPath, DUMPFOLDERPATH);
}

std::string getBackupPath(const std::string &name,
                          const std::string &dataPath) {
  return getFolderPath(name, dataPath, BACKUPFOLDERPATH);
}

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

// Helper function to copy a file to the backup folder if it's larger than a
// certain size
void createBackupOfFile(const fs::path &file, const fs::path &backupDir,
                        std::size_t minSize) {
  // Check the file size. No need to back up small files
  if (fs::file_size(file) > minSize) {
    fs::path destination = backupDir / file.filename();
    fs::create_directories(backupDir);
    // Check if there is a name conflict
    if (fs::exists(destination)) {
      // We now need to change the name so that the backup file is given a
      // unique name
      int counter = 1;
      std::string stem = destination.stem().string();
      std::string extension = destination.extension().string();
      // increment until a filename that doesn't exist is found
      while (true) {
        fs::path candidate =
            backupDir / (stem + "_" + std::to_string(counter) + extension);
        if (!fs::exists(candidate)) {
          destination = candidate;
          break;
        }
        ++counter;
      }
    }

    // Now copy the file to the final unique destination.
    fs::copy_file(file, destination, fs::copy_options::overwrite_existing);

    if (!isQuiet()) {
      std::cout << "Created backup of: " << file << std::endl;
    }
  } else {
    // File is too small to bother backing up
  }
}

// Function to back up a folder if it contains at least 10 files
void createBackupOfFolder(const fs::path &folder, const fs::path &backupDir) {
  // Check if the folder exists and is a directory
  if (!fs::exists(folder) || !fs::is_directory(folder)) {
    // std::cerr << "Folder does not exist or is not a directory: "
    //           << folder
    //           << std::endl;
    return;
  }

  // Count the number of files in the folder
  std::size_t fileCount = 0;
  for (const auto &entry : fs::directory_iterator(folder)) {
    if (fs::is_regular_file(entry)) {
      ++fileCount;
    }
  }

  // Proceed only if there are at least 10 files
  if (fileCount < 10) {
    // No need to warn about no backup if the folder is empty
    if (fileCount > 0) {
      if (!isQuiet()) {
        std::cout << "Folder contains fewer than 10 files. No backup created."
                  << std::endl;
      }
    }
    return;
  }

  // Prepare the destination folder path
  std::string folderName = folder.lexically_normal().filename();
  fs::path destination = backupDir / folderName;
  fs::create_directories(backupDir);

  // Resolve name conflict for the folder
  if (fs::exists(destination)) {
    int counter = 1;
    while (true) {
      fs::path candidate =
          backupDir / (folderName + "_" + std::to_string(counter));
      if (!fs::exists(candidate)) {
        destination = candidate;
        break;
      }
      ++counter;
    }
  }

  // Copy the entire folder to the backup directory
  try {
    fs::copy(folder, destination, fs::copy_options::recursive);
    if (!isQuiet()) {
      std::cout << "Created backup of folder: " << folder << " to "
                << destination << std::endl;
    }
  } catch (const fs::filesystem_error &e) {
    if (!isQuiet()) {
      std::cerr << "Error during folder backup: " << e.what() << std::endl;
    }
  }
}

// Clears a subfolder, copying large CSV files to a backup folder instead of
// deleting them
void clearOutputFolder(std::string name, std::string dataPath) {
  // Paths of folders where we delete all the files of specified extension, and
  // backup large csv files
  std::vector<std::string> paths = {getOutputPath(name, dataPath)};
  std::vector<std::string> extensionsToDelete = {".vtu", ".pvd", ".csv"};
  fs::path backupDir = getBackupPath(name, dataPath);

  for (const std::string &path : paths) {
    if (!fs::exists(path) || !fs::is_directory(path)) {
      continue;
    }

    for (const auto &entry : fs::directory_iterator(path)) {
      if (entry.is_regular_file()) {
        std::string extension = entry.path().extension().string();
        if (std::find(extensionsToDelete.begin(), extensionsToDelete.end(),
                      extension) != extensionsToDelete.end()) {
          if (extension == ".csv") {
            // Special handling for CSV files
            createBackupOfFile(entry.path(), backupDir, 100 * 1024); // 100KB
          }
          // Delete the file
          fs::remove(entry.path());
        }
      }
    }
  }

  // We also need to clear the data folder, but we always create a backup of
  // these if there are more than 10 vtu files
  createBackupOfFolder(getDataPath(name, dataPath), backupDir);
  fs::remove_all(getDataPath(name, dataPath));
}

// If we want to store some data that does not depend on either the node or
// cell, it is inefficient to store the data multiple times. The simplest way I
// have found to store extra data is by including it in the file name, dataPath.
// Example: The variable foo and bar are stored as "_foo=0.32_bar=4_".
std::string makeFileName(const Mesh &mesh, std::string name) {
  std::stringstream ss;
  ss << name << "_load=" << mesh.load
     << "_nrM=" << mesh.nr_elements_with_m3_fix_change << '_';
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
  std::vector<int> nodeIndexes(maxNrNodes);
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
} // namespace

std::string writeMeshToVtu(const Mesh &mesh, std::string folderName,
                           std::string dataPath, std::string fileName,
                           bool minimizationStep, VtuFieldLevel level,
                           std::string nameSuffix, std::string subDataFolder) {

  const int dim = 3;
  const int cell_size = 3;
  int nrElements = mesh.nrElements;

  // Due to periodic reconnecting, the number of nodes we use to represent the
  // periodic elements is not constant. We therefore create a node list here
  auto nodeMap = constructGhostNodeIndexes(mesh);
  int nrNodes = nodeMap.size();
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

  std::vector<double> force(nrNodes * dim);
  // boolean values represented by 0.0 and 1.0
  std::vector<double> fixed;
  std::vector<double> refIndex; // Index of reference node
  if (wantFixed) {
    fixed.resize(nrNodes);
  }
  if (wantRefIndex) {
    refIndex.resize(nrNodes);
  }
  std::vector<int> connectivity(nrElements * cell_size); // Mesh order
  std::vector<double> energy(nrElements);
  std::vector<double> det;
  std::vector<double> F11, F12, F21, F22;
  std::vector<double> F_Fix11, F_Fix12, F_Fix21, F_Fix22;
  std::vector<double> C11, C12, C22;
  std::vector<double> C_Fix11, C_Fix12, C_Fix22;
  std::vector<double> G11, G12, G22;
  std::vector<double> P11, P12, P21, P22;
  std::vector<double> sigma11, sigma12, sigma22;

  std::vector<double> largeAngle;
  std::vector<double> smallAngle;

  // These should be int, but the library i am using only takes doubles
  std::vector<double> nrm3;               // Int
  std::vector<double> deltaNrm3;          // Int
  std::vector<double> m11, m12, m21, m22; // Int
  std::vector<double> redQuadrant;        // Int
  std::vector<double> redQuadrantFixed;   // Int

  if (wantDet) {
    det.resize(nrElements);
  }
  if (wantMatrices) {
    F11.resize(nrElements);
    F12.resize(nrElements);
    F21.resize(nrElements);
    F22.resize(nrElements);
    F_Fix11.resize(nrElements);
    F_Fix12.resize(nrElements);
    F_Fix21.resize(nrElements);
    F_Fix22.resize(nrElements);
    C11.resize(nrElements);
    C12.resize(nrElements);
    C22.resize(nrElements);
    C_Fix11.resize(nrElements);
    C_Fix12.resize(nrElements);
    C_Fix22.resize(nrElements);
    G11.resize(nrElements);
    G12.resize(nrElements);
    G22.resize(nrElements);
    P11.resize(nrElements);
    P12.resize(nrElements);
    P21.resize(nrElements);
    P22.resize(nrElements);
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
  if (wantQuadrants) {
    redQuadrant.resize(nrElements);
    redQuadrantFixed.resize(nrElements);
  }

  leanvtk::VTUWriter writer;

  // Instead of getting the data directly from the nodes in the mesh, we extract
  // the data from the nodes in the elements in the mesh. This is because they
  // have a displaced position and to not result in overlapping elements
  // DO NOT USE std::vector<bool>! This leads to memory
  // coruption errors that are difficult to track down
  std::vector<char> alreadyCopied(nrNodes, false);
  // Iterate over each element in the mesh
  for (int elementIndex = 0; elementIndex < nrElements; ++elementIndex) {
    const TElement &e = mesh.elements[elementIndex];
    // Iterate over each node in the element
    for (size_t j = 0; j < e.ghostNodes.size(); ++j) {
      const GhostNode &gn = e.ghostNodes[j];
      // Element index
      auto it = nodeMap.find(gn.id);
      assert(it != nodeMap.end()); // or handle the error
      if (it == nodeMap.end()) {
        throw std::runtime_error("Ghost node not found in node map.");
      }
      int nodeIndex = it->second;
      if (!alreadyCopied[nodeIndex]) {
        points[nodeIndex * dim + 0] = gn.pos[0];
        points[nodeIndex * dim + 1] = gn.pos[1];
        points[nodeIndex * dim + 2] = 0;
        // We need to take the force from the mesh nodes, not the element nodes
        const Node &n = *mesh[gn];
        force[nodeIndex * dim + 0] = n.f[0];

        force[nodeIndex * dim + 1] = n.f[1];
        force[nodeIndex * dim + 2] = 0;
        if (wantFixed) {
          fixed[nodeIndex] = n.fixedNode;
        }
        if (wantRefIndex) {
          refIndex[nodeIndex] = gn.referenceId.i;
        }
        alreadyCopied[nodeIndex] = true;
      }
      // Here we define the meshing
      connectivity[elementIndex * cell_size + j] = nodeIndex;
    }

    energy[elementIndex] = e.energy;
    if (wantDet) {
      det[elementIndex] = e.C.determinant();
    }

    if (wantMatrices) {
      F11[elementIndex] = e.F(0, 0);
      F12[elementIndex] = e.F(0, 1);
      F21[elementIndex] = e.F(1, 0);
      F22[elementIndex] = e.F(1, 1);
      F_Fix11[elementIndex] = e.F_fixed_ref(0, 0);
      F_Fix12[elementIndex] = e.F_fixed_ref(0, 1);
      F_Fix21[elementIndex] = e.F_fixed_ref(1, 0);
      F_Fix22[elementIndex] = e.F_fixed_ref(1, 1);

      C11[elementIndex] = e.C(0, 0);
      C12[elementIndex] = e.C(0, 1);
      C22[elementIndex] = e.C(1, 1);
      C_Fix11[elementIndex] = e.C_fixed_ref(0, 0);
      C_Fix12[elementIndex] = e.C_fixed_ref(0, 1);
      C_Fix22[elementIndex] = e.C_fixed_ref(1, 1);

      G11[elementIndex] = e.G(0, 0);
      G12[elementIndex] = e.G(0, 1);
      G22[elementIndex] = e.G(1, 1);
      P11[elementIndex] = e.P(0, 0);
      P12[elementIndex] = e.P(0, 1);
      P21[elementIndex] = e.P(1, 0);
      P22[elementIndex] = e.P(1, 1);
      m11[elementIndex] = e.m(0, 0);
      m12[elementIndex] = e.m(0, 1);
      m21[elementIndex] = e.m(1, 0);
      m22[elementIndex] = e.m(1, 1);
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
    deltaNrm3[elementIndex] = e.m3Nr_fix - e.pastM3Nr_fix;
    if (wantQuadrants) {
      redQuadrant[elementIndex] = e.red_quadrant;
      redQuadrantFixed[elementIndex] = e.red_quadrant_fixed;
    }
  }

  // Debug confirm that all the nodes have been written to
  assert(std::all_of(alreadyCopied.begin(), alreadyCopied.end(),
                     [](bool value) { return value; }));

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
    writer.add_cell_scalar_field("F_Fix11", F_Fix11);
    writer.add_cell_scalar_field("F_Fix12", F_Fix12);
    writer.add_cell_scalar_field("F_Fix21", F_Fix21);
    writer.add_cell_scalar_field("F_Fix22", F_Fix22);
    writer.add_cell_scalar_field("C11", C11);
    writer.add_cell_scalar_field("C12", C12);
    writer.add_cell_scalar_field("C22", C22);
    writer.add_cell_scalar_field("C_Fix11", C_Fix11);
    writer.add_cell_scalar_field("C_Fix12", C_Fix12);
    writer.add_cell_scalar_field("C_Fix22", C_Fix22);
    writer.add_cell_scalar_field("G11", G11);
    writer.add_cell_scalar_field("G12", G12);
    writer.add_cell_scalar_field("G22", G22);
    writer.add_cell_scalar_field("P11", P11);
    writer.add_cell_scalar_field("P12", P12);
    writer.add_cell_scalar_field("P21", P21);
    writer.add_cell_scalar_field("P22", P22);
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
  if (wantQuadrants) {
    writer.add_cell_scalar_field("red_quadrant", redQuadrant);
    writer.add_cell_scalar_field("red_quadrant_fixed", redQuadrantFixed);
  }

  writer.add_vector_field("stress_field", force, dim);

  // write data
  writer.write_surface_mesh(filePath, dim, cell_size, points, connectivity);
  // Save a compressed version
  // Paraview doesn't open compressed files >:(
  // TODO, if used, remember to update the dataPath properly so pvd collection
  // is updated.
  // filePath = compressFile(filePath);
  return filePath;
}

// Function to join strings with a delimiter
std::string join(const std::vector<std::string> &strings,
                 const std::string &delimiter) {
  std::string result;
  for (size_t i = 0; i < strings.size(); ++i) {
    result += strings[i];
    if (i < strings.size() - 1)
      result += delimiter;
  }
  return result;
}

// Function to initialize a CSV file for writing
std::ofstream initCsvFile(const std::string &folderName,
                          const std::string &dataPath, const Simulation &s,
                          const std::string subFolder) {
  // Construct the full file path
  std::string fileName = MACRODATANAME + std::string(".csv");
  std::string outputPath = getFolderPath(folderName, dataPath, subFolder);
  fs::path filePath = fs::path(outputPath) / fileName;

  bool headerWasWritten = insertHeaderIfNeeded(filePath, s);

  if (!headerWasWritten) {
    // If we start from a dump, we need to trim the csv file
    // Before we do, we will create a backup of the file
    fs::path backupDir = getBackupPath(s.simName, dataPath);
    // create backup if it is larger than 100KB
    createBackupOfFile(filePath, backupDir, 100 * 1024); // 100KB
    trimCsvFile(filePath, s);
  }

  // Open the file in append mode
  std::ofstream file(filePath, std::ios::app);
  if (!file.is_open()) {
    // Handle the error if file cannot be opened
    throw std::runtime_error("Unable to open file: " + filePath.string());
  }

  // Construct file path
  return file;
}

// Helper function to split a line by commas
std::vector<std::string> splitLine(const std::string &line) {
  std::stringstream ss(line);
  std::string item;
  std::vector<std::string> elements;
  while (std::getline(ss, item, ',')) {
    elements.push_back(item);
  }
  return elements;
}

struct CsvRemap {
  std::vector<int> srcIndexForDst;
  size_t loadStepIndex = std::string::npos;
};

CsvRemap buildCsvRemap(const std::vector<std::string> &fileHeaders,
                       const std::vector<std::string> &expectedHeaders) {
  CsvRemap remap;
  remap.srcIndexForDst.assign(expectedHeaders.size(), -1);
  for (size_t i = 0; i < expectedHeaders.size(); ++i) {
    const auto &name = expectedHeaders[i];
    if (name == "load_step") {
      remap.loadStepIndex = i;
    }
    for (size_t j = 0; j < fileHeaders.size(); ++j) {
      if (fileHeaders[j] == name) {
        remap.srcIndexForDst[i] = static_cast<int>(j);
        break;
      }
    }
  }
  return remap;
}

std::vector<std::string> remapCsvRow(const std::vector<std::string> &row,
                                     const CsvRemap &remap,
                                     size_t expectedSize) {
  std::vector<std::string> out(expectedSize, "");
  for (size_t i = 0; i < expectedSize; ++i) {
    int srcIndex = remap.srcIndexForDst[i];
    if (srcIndex >= 0 && static_cast<size_t>(srcIndex) < row.size()) {
      out[i] = row[static_cast<size_t>(srcIndex)];
    }
  }
  return out;
}

/*
  This function trims any rows from the CSV file whose loadStep
  is larger than the current simulation's loadStep. It stops at
  the first row that meets that condition and removes all
  subsequent rows (including that one).
*/
void trimCsvFile(const std::string &filePath, const Simulation &s) {
  // The current loadStep from the simulation
  int currentLoadStep = s.mesh.loadSteps;

  // Open the CSV file for reading
  std::ifstream inputFile(filePath);
  if (!inputFile.is_open()) {
    throw std::runtime_error("Could not open file for reading: " + filePath);
  }

  // We'll store the lines we want to keep in this vector
  std::vector<std::string> lines;
  std::string line;

  // This flag will tell us if we found a line with a loadStep
  // larger (or equal) to the current loadStep
  bool foundLargerStep = false;
  bool firstLine = true;
  std::vector<std::string> expectedHeaders = getCsvHeaders(s);
  CsvRemap remap;

  // Read line by line
  while (std::getline(inputFile, line)) {
    // Split the line into columns
    std::vector<std::string> elements = splitLine(line);

    if (firstLine) {
      remap = buildCsvRemap(elements, expectedHeaders);
      if (remap.loadStepIndex == std::string::npos ||
          remap.srcIndexForDst[remap.loadStepIndex] < 0) {
        throw std::runtime_error(
            "CSV header missing required load_step column: " + filePath);
      }
      lines.push_back(join(expectedHeaders, ","));
      firstLine = false;
      continue;
    }

    // Not first line
    std::vector<std::string> remapped =
        remapCsvRow(elements, remap, expectedHeaders.size());

    try {
      // Convert the relevant column to an integer
      const std::string &loadStepStr = remapped[remap.loadStepIndex];
      if (loadStepStr.empty()) {
        throw std::runtime_error("Missing load_step value in line: " + line);
      }
      int lineLoadStep = std::stoi(loadStepStr);

      // If the line's loadStep is >= the current loadStep,
      // we stop reading further and break out of the loop
      if (lineLoadStep > currentLoadStep) {
        foundLargerStep = true;
        break;
      }
    } catch (const std::invalid_argument &e) {
      // If we can't parse an integer, decide how you want to handle it
      throw std::runtime_error("Invalid loadStep value in line: " + line);
    } catch (const std::out_of_range &e) {
      throw std::runtime_error("LoadStep value out of range in line: " + line);
    }

    // If we haven't found a larger step, keep the line
    lines.push_back(join(remapped, ","));
  }

  // Close the input file
  inputFile.close();

  // If we never found a line with loadStep >= current loadStep,
  // then there's nothing to trim, so we simply return.
  if (!foundLargerStep) {
    return;
  }

  // Otherwise, rewrite the file with only the kept lines
  std::ofstream outputFile(filePath, std::ios::trunc);
  if (!outputFile.is_open()) {
    throw std::runtime_error("Could not open file for writing: " + filePath);
  }

  // Write back the lines we decided to keep
  for (const auto &l : lines) {
    outputFile << l << "\n";
  }

  // Automatically closes when going out of scope
}

void saveConfigFile(Config conf, std::string dataPath) {
  // Construct the full file path
  fs::path filePath = fs::path(getOutputPath(conf.name, dataPath)) / CONFIGNAME;

  // Check if both paths exist before checking equivalence
  if (fs::exists(conf.configPath) && fs::exists(filePath)) {
    // Check if the source and destination paths are the same
    if (fs::equivalent(conf.configPath, filePath)) {
      return;
    }
  }

  try {
    std::ifstream src(conf.configPath, std::ios::binary);
    if (!src) {
      if (!isQuiet()) {
        if (conf.configPath.empty()) {
          std::cout << "No config path specified." << std::endl;
        } else {
          std::cerr << "Failed to open source file: " << conf.configPath
                    << std::endl;
        }
      }
      return;
    }

    // Prepare the destination
    std::ofstream dst(filePath, std::ios::binary);

    if (!dst) {
      if (!isQuiet()) {
        std::cerr << "Failed to open destination file: " << filePath
                  << std::endl;
        std::cerr << "Check if the output directory exists and you have "
                     "permission to write."
                  << std::endl;
      }
      return;
    }

    dst << src.rdbuf(); // Copy the content
  } catch (const std::exception &e) {
    if (!isQuiet()) {
      std::cerr << "Exception occurred: " << e.what() << std::endl;
    }
  }
}

void saveConfigFile(std::string configFile, std::string dataPath) {
  Config conf = parseConfigFile(configFile);
  saveConfigFile(conf, dataPath);
}

// Function to write line to CSV using an open file stream
void writeLineToCsv(std::ofstream &file,
                    const std::vector<std::string> &strings) {
  if (!file.is_open()) {
    throw std::runtime_error("File stream is not open.");
  }

  // Join the strings into a single line
  std::string line = join(strings, ",");

  // Write the line to file
  file << line << '\n';

  // Check if write was successful
  if (!file) {
    throw std::runtime_error("Failed to write to file.");
  }
}

void writeLineToCsv(std::ofstream &file, const std::vector<double> &values) {
  std::vector<std::string> stringValues;
  stringValues.reserve(values.size());
  for (double value : values) {
    stringValues.push_back(std::to_string(value));
  }
  writeLineToCsv(file, stringValues);
}

// This writes any information we want to one line of the csv file
void writeToCsv(std::ofstream &file, const Simulation &s) {
  const auto lineData = getStringVector(s);
  writeLineToCsv(file, lineData);
}

void writeCsvHeaders(std::ofstream &file, const Simulation &s) {
  const auto lineData = getCsvHeaders(s);
  writeLineToCsv(file, lineData);
}

std::vector<std::string> getCsvHeaders(const Simulation &s) {
  const auto &cols = s.getCsvColumns();
  std::vector<std::string> headers;
  headers.reserve(cols.size());
  for (const auto &col : cols) {
    headers.push_back(col.name);
  }
  return headers;
}

std::vector<std::string> getStringVector(const Simulation &s) {
  const auto &cols = s.getCsvColumns();
  std::vector<std::string> row;
  row.reserve(cols.size());
  for (const auto &col : cols) {
    row.push_back(col.getter(s));
  }
  return row;
}

bool insertHeaderIfNeeded(const std::string &filename, const Simulation &s) {
  const fs::path p = fs::path(filename);

  // If file exists and has content, do nothing.
  if (fs::exists(p)) {
    // If the file exists but is empty, we still want to write the header.
    std::error_code ec;
    const auto sz = fs::file_size(p, ec);
    if (!ec && sz > 0) {
      return false;
    }
  }

  // Create (or truncate an empty) file and write header.
  std::ofstream fileOut(filename, std::ios::trunc);
  if (!fileOut.is_open()) {
    throw std::runtime_error("Unable to create/open file with header: " +
                             filename);
  }

  writeCsvHeaders(fileOut, s);
  return true;
}

void createCollection(const std::string &folderPath,
                      const std::string &destination,
                      const std::string &regexPattern,
                      const std::string &extension,
                      const std::vector<double> &timestep) {

  // Convert input paths to absolute to avoid ambiguities.
  fs::path absFolder = fs::absolute(folderPath);
  fs::path absDestination = fs::absolute(destination);

  // Check that the input folder exists.
  if (!fs::exists(absFolder)) {
    throw std::runtime_error("The folder path does not exist: " +
                             absFolder.string());
  }

  // Container for storing files along with the extracted number from their
  // filenames.
  std::vector<std::pair<int, fs::path>> filesWithNumbers;
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
      std::smatch match;
      if (std::regex_match(filename, match, regex) && match.size() == 2) {
        int number = std::stoi(match[1].str());
        filesWithNumbers.emplace_back(number, entry.path());
      }
    }
  }

  // Sort the files based on the extracted numeric value.
  std::sort(filesWithNumbers.begin(), filesWithNumbers.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  // Construct the full path of the collection file (output file).
  fs::path collectionFilePath =
      absDestination / (std::string(COLLECTIONNAME) + ".pvd");
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
  for (std::size_t i = 0; i < filesWithNumbers.size(); ++i) {
    // Use the provided timestep if available, otherwise use the extracted file
    // number.
    double timeValue =
        (i < timestep.size()) ? timestep[i] : filesWithNumbers[i].first;

    // Compute the absolute file path and then its relative path from the
    // destination directory.
    fs::path absFilePath = fs::absolute(filesWithNumbers[i].second);
    fs::path fileRelativePath = fs::relative(absFilePath, absDestination);

    // Write the DataSet element, ensuring the file attribute is quoted.
    outFile << "    <DataSet timestep=\"" << timeValue
            << "\" group=\"\" part=\"0\" file=\"" << fileRelativePath.string()
            << "\"/>\n";
  }

  // Close the XML tags.
  outFile << "  </Collection>\n";
  outFile << "</VTKFile>\n";

  outFile.close();
}

bool simulationAlreadyComplete(const std::string &name,
                               const std::string &dataPath, double maxLoad) {
  std::string filePath = getOutputPath(name, dataPath) + MACRODATANAME + ".csv";

  std::ifstream file(filePath);
  if (!file) {
    return false; // File couldn't be opened
  }
  std::string line, header;
  if (!std::getline(file, header)) {
    return false; // Empty file
  }

  // Find "load" column index
  std::istringstream headerStream(header);
  std::string column;
  size_t loadIndex = std::string::npos, index = 0;

  while (std::getline(headerStream, column, ',')) {
    if (column == "load") {
      loadIndex = index;
      break;
    }
    index++;
  }

  if (loadIndex == std::string::npos) {
    if (!isQuiet()) {
      std::cerr << "Error! Load not found in macroData file!\n";
    }
    return false; // "Load" column missing
  }

  // Read last non-empty line
  std::string lastLine;
  while (std::getline(file, line)) {
    if (!line.empty()) {
      lastLine = line;
    }
  }

  if (lastLine.empty()) {
    return false; // No data rows
  }

  // Extract last row's load value
  std::istringstream lastLineStream(lastLine);
  index = 0;
  double lastLoad = 0.0;

  try {
    while (std::getline(lastLineStream, column, ',')) {
      if (index == loadIndex) {
        lastLoad = std::stod(column);
        break;
      }
      index++;
    }
  } catch (const std::exception &) {
    return false; // Handle conversion failure
  }

  return lastLoad >= maxLoad;
}
