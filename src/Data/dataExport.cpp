#include "dataExport.h"

#include "../Simulation/simulation.h"
#include "Data/lean_vtk.h"
#include "Data/paramParser.h"
#include "settings.h"
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <iostream>
#include <iterator>

std::string findOutputPath() {
  // Define the paths to check
  std::vector<std::string> paths = {"/Volumes/data/",
                                    "/media/elias/dataStorage/",
                                    "/data2/elundheim/", "/data/elundheim/"};

  // Initialize a variable to store the chosen path
  std::string chosen_path;

  // Iterate through the paths and check if they exist
  for (const auto &path : paths) {
    if (std::filesystem::exists(path)) {
      chosen_path = path;
      break; // Stop the loop once a valid path is found
    }
  }

  // Check if a valid path was found or throw an error
  if (chosen_path.empty()) {
    throw std::runtime_error("None of the provided paths exist.");
  } else {
    // We now also add the output folder name
    chosen_path += OUTPUTFOLDERPATH;
    std::cout << "Chosen output path: " << chosen_path << std::endl;
  }

  return chosen_path;
}

std::string searchForConfig(std::string dumpPath) {
  // Extract the dataFolder path
  std::string dataFolderPath =
      dumpPath.substr(0, dumpPath.rfind(DUMPFOLDERPATH));

  // Append config.conf to the dataFolder path
  std::string configFilePath = dataFolderPath + CONFIGNAME;

  // Check if the config.conf file exists
  if (std::filesystem::exists(configFilePath)) {
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

namespace fs = std::filesystem;

std::string getFolderPath(const std::string &name, const std::string &dataPath,
                          const fs::path &subfolder = "") {
  return (fs::path(dataPath) / name / subfolder).string() + '/';
}

std::string getOutputPath(const std::string &name,
                          const std::string &dataPath) {
  return getFolderPath(name, dataPath);
}

std::string getDataPath(const std::string &name, const std::string &dataPath) {
  return getFolderPath(name, dataPath, DATAFOLDERPATH);
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

void createDataFolder(std::string name, std::string dataPath) {
  std::vector<std::string> paths = {getDataPath(name, dataPath),
                                    getFramePath(name, dataPath),
                                    getDumpPath(name, dataPath)};
  for (std::string path : paths) {
    // Ensure the directory exists
    std::filesystem::create_directories(path);
  }
}

// Get file path, and check that the path exsists
std::string getFilePath(std::string fileName, std::string folderName,
                        std::string dataPath, std::string fileType = ".vtu") {
  std::string fileNameWithType = fileName + fileType;
  std::string directory = getDataPath(folderName, dataPath);

  // Check if the directory exists
  if (!std::filesystem::exists(directory) ||
      !std::filesystem::is_directory(directory)) {
    throw std::runtime_error("Directory does not exist: " + directory +
                             "\nHave you run the createDataFolder function?");
  }

  return directory + fileNameWithType;
}

// Helper function to copy a file to the backup folder if it's larger than a
// certain size
void backupLargeFile(const std::filesystem::path &file,
                     const std::filesystem::path &backupDir,
                     std::size_t maxSize) {
  // Check the file size
  if (std::filesystem::file_size(file) > maxSize) {
    std::cout << "large file\n";
    // Ensure the backup directory exists
    std::filesystem::create_directories(backupDir);
    // Construct the destination path
    std::filesystem::path destination = backupDir / file.filename();
    // Copy the file
    std::cout << backupDir << destination;
    std::filesystem::copy_file(
        file, destination, std::filesystem::copy_options::overwrite_existing);
  }
}

// Clears a subfolder, copying large CSV files to a backup folder instead of
// deleting them
void clearOutputFolder(std::string name, std::string dataPath) {
  std::vector<std::string> paths = {getOutputPath(name, dataPath),
                                    getDataPath(name, dataPath),
                                    getFramePath(name, dataPath)};
  std::vector<std::string> extensionsToDelete = {".vtu", ".pvd", ".csv", ".png",
                                                 ".log"};
  std::filesystem::path backupDir = getBackupPath(name, dataPath);

  for (const std::string &path : paths) {
    if (!std::filesystem::exists(path) ||
        !std::filesystem::is_directory(path)) {
      continue;
    }

    for (const auto &entry : std::filesystem::directory_iterator(path)) {
      if (entry.is_regular_file()) {
        std::string extension = entry.path().extension().string();
        if (std::find(extensionsToDelete.begin(), extensionsToDelete.end(),
                      extension) != extensionsToDelete.end()) {
          if (extension == ".csv") {
            // Special handling for CSV files
            backupLargeFile(entry.path(), backupDir, 100 * 1024); // 100KB
          }
          // Delete the file
          std::filesystem::remove(entry.path());
        }
      }
    }
  }
}

// If we want to store some data that does not depend on either the node or
// cell, it is inefficient to store the data multiple times. The simplest way I
// have found to store extra data is by including it in the file name, dataPath.
// Example: The variable foo and bar are stored as "_foo=0.32_bar=4_".
std::string makeFileName(const Mesh &mesh, std::string name) {
  std::stringstream ss;
  ss << name << "_load=" << mesh.load << "_nrM=" << mesh.nrPlasticChanges
     << '_';
  return ss.str();
}

void writeMeshToVtu(const Mesh &mesh, std::string folderName,
                    std::string dataPath) {

  const int dim = 3;
  const int cell_size = 3;
  int n = mesh.rows;
  int m = mesh.cols;
  if (mesh.usingPBC) {
    // We create an extra row and column for the ghost nodes.
    n += 1;
    m += 1;
  }
  int nm = n * m;
  // Since timeStep is static, it will increase each time the function is
  // called.
  int nrNodes = nm;
  int nrElements = mesh.nrElements;

  std::string fileName = makeFileName(mesh, folderName);

  std::string filePath;

  filePath = getFilePath(fileName + "." + std::to_string(mesh.loadSteps),
                         folderName, dataPath);

  std::vector<double> points(nrNodes * dim);
  std::vector<double> force(nrNodes * dim);
  std::vector<double> fixed(
      nrNodes); // boolean values represented by 0.0 and 1.0
  std::vector<int> elements(nrElements * cell_size);
  std::vector<double> energy(nrElements);
  std::vector<double> C11(nrElements);
  std::vector<double> C12(nrElements);
  std::vector<double> C22(nrElements);
  std::vector<double> P11(nrElements);
  std::vector<double> P12(nrElements);
  std::vector<double> P21(nrElements);
  std::vector<double> P22(nrElements);
  std::vector<double> resolvedShearStress(nrElements);

  leanvtk::VTUWriter writer;

  // Instead of getting the data directly from the nodes in the mesh, we extract
  // the data from the nodes in the elements in the mesh. This is because they
  // have a displaced position and to not result in overlapping elements
  std::vector<char> alreadyCopied(
      nrNodes, false); // DO NOT USE std::vector<bool>! This leads to memory
                       // coruption errors that are difficult to track down
  // Iterate over each element in the mesh
  for (int elementIndex = 0; elementIndex < nrElements; ++elementIndex) {
    const TElement &e = mesh.elements[elementIndex];
    // Iterate over each node in the element
    for (size_t j = 0; j < e.nodes.size(); ++j) {
      const Node &n = e.nodes[j];
      // Element index
      int nodeIndex = mesh.usingPBC ? n.ghostId.i : n.id.i;
      if (!alreadyCopied[nodeIndex]) {
        points[nodeIndex * dim + 0] = n.pos()[0];
        points[nodeIndex * dim + 1] = n.pos()[1];
        points[nodeIndex * dim + 2] = 0;
        force[nodeIndex * dim + 0] = n.f[0];
        force[nodeIndex * dim + 1] = n.f[1];
        force[nodeIndex * dim + 2] = 0;
        fixed[nodeIndex] = n.fixedNode;
        alreadyCopied[nodeIndex] = true;
      }

      // We choose to either use the ghost id or the real id depending on
      // whether or not we are using pbc.
      elements[elementIndex * cell_size + j] = nodeIndex;
    }

    energy[elementIndex] = e.energy;
    C11[elementIndex] = e.C(0, 0);
    C12[elementIndex] = e.C(0, 1);
    C22[elementIndex] = e.C(1, 1);
    P11[elementIndex] = e.P(0, 0);
    P12[elementIndex] = e.P(0, 1);
    P21[elementIndex] = e.P(1, 0);
    P22[elementIndex] = e.P(1, 1);
    resolvedShearStress[elementIndex] = e.resolvedShearStress;
  }

  // Debug confirm that all the nodes have been written to
  assert(std::all_of(alreadyCopied.begin(), alreadyCopied.end(),
                     [](bool value) { return value; }));

  // connect data to writer
  writer.add_cell_scalar_field("energy_field", energy);
  writer.add_cell_scalar_field("resolvedShearStress", resolvedShearStress);
  writer.add_scalar_field("fixed", fixed);
  writer.add_cell_scalar_field("C11", C11);
  writer.add_cell_scalar_field("C12", C12);
  writer.add_cell_scalar_field("C22", C22);
  writer.add_cell_scalar_field("P11", P11);
  writer.add_cell_scalar_field("P12", P12);
  writer.add_cell_scalar_field("P21", P21);
  writer.add_cell_scalar_field("P22", P22);

  writer.add_vector_field("stress_field", force, dim);

  // write data
  writer.write_surface_mesh(filePath, dim, cell_size, points, elements);
}

#include <fstream>
#include <string>
#include <vector>

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
                          const std::string &dataPath) {
  // Construct the full file path
  std::string filePath =
      getOutputPath(folderName, dataPath) + MACRODATANAME + ".csv";

  insertHeaderIfNeeded(filePath);

  // Open the file in append mode
  std::ofstream file(filePath, std::ios::app);
  if (!file.is_open()) {
    // Handle the error if file cannot be opened
    throw std::runtime_error("Unable to open file: " + filePath);
  }

  // Construct file path
  return file;
}

// Duplicate the config file given to the new output folder
void saveConfigFile(Config conf) {
  // Construct the full file path
  std::string filePath =
      getOutputPath(conf.name, findOutputPath()) + CONFIGNAME;

  std::ifstream src(conf.configPath,
                    std::ios::binary); // Open the source file in binary mode
  std::ofstream dst(
      filePath, std::ios::binary); // Open the destination file in binary mode

  if (!src) {
    std::cerr << "Failed to open source file: " << conf.configPath << std::endl;
    return;
  }

  if (!dst) {
    std::cerr << "Failed to open destination file: " << filePath << std::endl;
    std::cerr << "Did you remember to first create the folder?" << std::endl;
    return;
  }

  dst << src.rdbuf(); // Copy the content
  src.close();
  dst.close();
}

void saveConfigFile(std::string configFile) {
  Config conf = parseConfigFile(configFile);
  saveConfigFile(conf);
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
  file << line << std::endl;
}

void writeLineToCsv(std::ofstream &file, const std::vector<double> &values) {
  std::vector<std::string> stringValues;
  stringValues.reserve(values.size());
  for (double value : values) {
    stringValues.push_back(std::to_string(value));
  }
  writeLineToCsv(file, stringValues);
}

// ChatGPT magic. Converts everything into strings
template <typename... Args>
std::vector<std::string> createStringVector(Args &&...args) {
  std::vector<std::string> vec;
  (vec.push_back([=] {
    std::ostringstream oss;
    oss << args;
    return oss.str();
  }()),
   ...);
  return vec;
}
// This writes any information we want to one line of the cvs file
void writeToCsv(std::ofstream &file, const Simulation &s) {
  static int lineCount = 0;
  lineCount += 1;
  // Must be in the same order as getCsvCols
  auto lineData = createStringVector(
      lineCount, s.mesh.load, s.mesh.averageEnergy, s.mesh.maxEnergy,
      s.mesh.averageResolvedShearStress(), s.mesh.nrPlasticChanges,
      s.FIRERep.nrIter, s.LBFGSRep.nrIter, s.FIRERep.nfev, s.LBFGSRep.nfev,
      s.FIRERep.terminationType, s.LBFGSRep.terminationType, s.getRunTime(),
      s.getEstimatedRemainingTime(), s.mesh.bounds[0], s.mesh.bounds[1],
      s.mesh.bounds[2], s.mesh.bounds[3], s.config.dtStart);

  writeLineToCsv(file, lineData);
}

std::vector<std::string> getCsvCols() {
  return {"Line nr",
          "Load",
          "Avg energy",
          "Max energy",
          "Avg RSS",
          "Nr plastic deformations",
          "Nr FIRE iterations",
          "Nr LBFGS iterations",
          "Nr FIRE func evals",
          "Nr LBFGS func evals",
          "FIRE Term reason",
          "LBFGS Term reason",
          "Run time",
          "Est time remaining",
          "maxX",
          "minX",
          "maxY",
          "minY",
          "dt_start"};
}

void writeCsvCols(std::ofstream &file) {
  auto lineData = getCsvCols();
  writeLineToCsv(file, lineData);
}

// Insert header into a CSV file if needed
void insertHeaderIfNeeded(const std::string &filename) {
  std::ifstream fileIn(filename);

  if (!fileIn.is_open()) {
    // Create the file with only the header if it does not exist
    std::ofstream fileOut(filename);
    if (!fileOut.is_open()) {
      throw std::runtime_error("Unable to create file.");
    }
    writeCsvCols(fileOut);
    fileOut.close();
    return;
  }

  // Read the first line from the file
  std::string firstLine;
  std::getline(fileIn, firstLine);
  fileIn.close();
  if (firstLine.empty()) {
    std::ofstream fileOut(filename);
    writeCsvCols(fileOut);
    fileOut.close();
  }
  // Check if the header needs to be updated
  // Here we do a super lazy check. We just assume that if the first line starts
  // with the same character as the first character of the first column, things
  // are as they should be.
  else if (firstLine[0] != getCsvCols()[0][0]) {
    // Create a temporary file for safe writing
    std::string tempFilename = filename + ".tmp";
    std::ofstream fileOut(tempFilename);
    if (!fileOut.is_open()) {
      throw std::runtime_error("Unable to open temporary file for writing.");
    }

    // Write the correct header
    writeCsvCols(fileOut);

    fileIn.open(filename); // Re-open original file to copy remaining lines
    std::string line;
    while (std::getline(fileIn, line)) {
      fileOut << line << std::endl;
    }
    fileIn.close();
    fileOut.close();

    // Replace the original file with the new one
    std::remove(filename.c_str());
    std::rename(tempFilename.c_str(), filename.c_str());
  }
}

void createCollection(const std::string folderPath,
                      const std::string destination,
                      const std::string collectionName,
                      const std::string extension,
                      const std::vector<double> &timestep) {
  using namespace std::filesystem;

  std::vector<std::pair<int, path>> filesWithNumbers;

  std::regex regexPattern(".*\\.([0-9]+)\\.vtu");

  for (const auto &entry : directory_iterator(folderPath)) {
    if (entry.path().extension() == extension) {
      std::smatch match;
      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, match, regexPattern) &&
          match.size() == 2) {
        int number = std::stoi(match[1].str());
        filesWithNumbers.emplace_back(number, entry.path());
      }
    }
  }

  // Sort the files based on the extracted number
  std::sort(filesWithNumbers.begin(), filesWithNumbers.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

  path relativePath = relative(
      folderPath,
      path(folderPath).parent_path()); // Get the relative path from the parent

  std::ofstream outFile(destination + "/" + collectionName + ".pvd");
  outFile << "<?xml version=\"1.0\"?>\n";
  outFile << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
  outFile << "<Collection>\n";

  for (size_t i = 0; i < filesWithNumbers.size(); ++i) {
    double ts = timestep.size() > i ? timestep[i] : static_cast<double>(i);
    outFile << "<DataSet timestep=\"" << ts
            << "\" group=\"\" part=\"0\" file=\"" << folderPath
            << filesWithNumbers[i].second.filename().string() << "\"/>\n";
  }

  outFile << "</Collection>\n";
  outFile << "</VTKFile>\n";
  outFile.close();
}

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

bool simulationAlreadyComplete(std::string name, std::string dataPath,
                               double maxLoad) {
  // Construct the full file path
  std::string filePath = getOutputPath(name, dataPath) + MACRODATANAME + ".csv";

  std::ifstream file(filePath);
  if (!file.is_open()) {
    // Handle file not found or cannot be opened
    return false;
  }

  std::string line;
  std::string header;
  std::string lastLine;

  // Read the first line (header)
  if (!std::getline(file, header)) {
    // Handle empty file or error reading header
    file.close();
    return false;
  }

  std::vector<std::string> headerColumns;
  std::stringstream headerStream(header);
  std::string column;
  while (std::getline(headerStream, column, ',')) {
    headerColumns.push_back(column);
  }

  auto loadIt = std::find(headerColumns.begin(), headerColumns.end(), "Load");
  if (loadIt == headerColumns.end()) {
    // Handle case where "Load" column is not found
    file.close();
    return false;
  }

  size_t loadIndex = std::distance(headerColumns.begin(), loadIt);

  // Read until the end to get the last line
  while (std::getline(file, line)) {
    if (!line.empty()) {
      lastLine = line;
    }
  }

  file.close();

  if (lastLine.empty()) {
    // If the file does not contain any data rows
    return false;
  }

  // Parse the last row
  std::stringstream lastLineStream(lastLine);
  std::vector<std::string> lastLineCells;
  while (std::getline(lastLineStream, column, ',')) {
    lastLineCells.push_back(column);
  }

  if (lastLineCells.size() <= loadIndex) {
    // If the last row does not have enough columns
    return false;
  }

  double lastLoad = std::stod(lastLineCells[loadIndex]);

  return lastLoad == maxLoad;
}