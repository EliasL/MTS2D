#include "dataExport.h"

#include "../Simulation/simulation.h"
#include "Data/lean_vtk.h"
#include "Data/paramParser.h"
#include "settings.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <ostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem; // Alias for filesystem

std::string findOutputPath() {
  // Define the paths to check
  std::vector<std::string> paths = {
      "/Volumes/data/",
      "/media/elias/dataStorage/",
      "/data2/elundheim/", // PMMH
      "/data/elundheim/",  // PMMH
      "/Users/elias/Work/PhD/Code/localData/",
      "/lustre/fswork/projects/rech/bph/uog82gz/", // JeanZay
  };

  // Initialize a variable to store the chosen path
  std::string chosen_path;

  // Iterate through the paths and check if they exist
  for (const auto &path : paths) {
    if (fs::exists(path)) {
      chosen_path = path;
      break; // Stop the loop once a valid path is found
    }
  }

  // Check if a valid path was found or throw an error
  if (chosen_path.empty()) {
    throw std::runtime_error(
        "Out path does not exsist. Is your storage device connected?");
  } else {
    // We now also add the output folder name
    chosen_path += OUTPUTFOLDERPATH;
    // std::cout << "Chosen output path: " << chosen_path << std::endl;
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

namespace fs = fs;

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
    fs::create_directories(path);
  }
}

// Get file path, and check that the path exsists
std::string getFilePath(std::string fileName, std::string folderName,
                        std::string dataPath, std::string fileType = ".vtu") {
  std::string fileNameWithType = fileName + fileType;
  std::string directory = getDataPath(folderName, dataPath);

  // Check if the directory exists
  if (!fs::exists(directory) || !fs::is_directory(directory)) {
    throw std::runtime_error("Directory does not exist: " + directory +
                             "\nHave you run the createDataFolder function?");
  }

  return directory + fileNameWithType;
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

    std::cout << "Created backup of: " << file << std::endl;
  } else {
    // File is too small to bother backing up
  }
}

// Function to back up a folder if it contains at least 10 files
void createBackupOfFolder(const fs::path &folder, const fs::path &backupDir) {
  // Check if the folder exists and is a directory
  if (!fs::exists(folder) || !fs::is_directory(folder)) {
    std::cerr << "Error: Folder does not exist or is not a directory."
              << std::endl;
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
    std::cout << "Folder contains fewer than 10 files. No backup created."
              << std::endl;
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
    std::cout << "Created backup of folder: " << folder << " to " << destination
              << std::endl;
  } catch (const fs::filesystem_error &e) {
    std::cerr << "Error during folder backup: " << e.what() << std::endl;
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
  // boolean values represented by 0.0 and 1.0
  std::vector<double> fixed(nrNodes);
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
  std::vector<double> nrm1(nrElements);
  std::vector<double> nrm2(nrElements);
  std::vector<double> nrm3(nrElements);
  std::vector<double> m11(nrElements);
  std::vector<double> m12(nrElements);
  std::vector<double> m21(nrElements);
  std::vector<double> m22(nrElements);

  leanvtk::VTUWriter writer;

  // Instead of getting the data directly from the nodes in the mesh, we extract
  // the data from the nodes in the elements in the mesh. This is because they
  // have a displaced position and to not result in overlapping elements
  std::vector<char> alreadyCopied(
      nrNodes, false); // DO NOT USE std::vector<bool>! This leads to memory
                       // coruption errors that are difficult to track down
                       // (with my skills at least)
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
    nrm1[elementIndex] = e.m1Nr;
    nrm2[elementIndex] = e.m2Nr;
    nrm3[elementIndex] = e.m3Nr;
    m11[elementIndex] = e.m(0, 0);
    m12[elementIndex] = e.m(0, 1);
    m21[elementIndex] = e.m(1, 0);
    m22[elementIndex] = e.m(1, 1);
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

  writer.add_cell_scalar_field("nrm1", nrm1);
  writer.add_cell_scalar_field("nrm2", nrm2);
  writer.add_cell_scalar_field("nrm3", nrm3);
  writer.add_cell_scalar_field("m11", m11);
  writer.add_cell_scalar_field("m12", m12);
  writer.add_cell_scalar_field("m21", m21);
  writer.add_cell_scalar_field("m22", m22);

  writer.add_vector_field("stress_field", force, dim);

  // write data
  writer.write_surface_mesh(filePath, dim, cell_size, points, elements);
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
                          const std::string &dataPath, const Simulation &s) {
  // Construct the full file path
  std::string filePath =
      getOutputPath(folderName, dataPath) + MACRODATANAME + ".csv";

  bool headerWasWritten = insertHeaderIfNeeded(filePath);

  if (!headerWasWritten) {
    // If we start from a dump, we need to trim the csv file
    // Before we do, we will create a backup of the file
    fs::path backupDir = getBackupPath(s.name, dataPath);
    createBackupOfFile(filePath, backupDir, 100 * 1024); // 100KB
    trimCsvFile(filePath, s);
  }

  // Open the file in append mode
  std::ofstream file(filePath, std::ios::app);
  if (!file.is_open()) {
    // Handle the error if file cannot be opened
    throw std::runtime_error("Unable to open file: " + filePath);
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

// When a simulation is resumed from a dump, unless the program was stopped
// right after the dump was created, the csv file will have lines that need to
// be overwritten
void trimCsvFile(const std::string &filePath, const Simulation &s) {
  int loadIndex = 0; // Index for load value in CSV
  double currentLoad = s.mesh.load;

  std::ifstream inputFile(filePath); // Open the CSV file for reading
  if (!inputFile.is_open()) {
    throw std::runtime_error("Could not open file for reading");
  }

  std::vector<std::string> lines; // Temporary storage for lines to keep
  std::string line;
  bool foundLargerLoad = false;

  // Read the file line by line
  while (std::getline(inputFile, line)) {
    lines.push_back(line); // Add line to the buffer if no match is found yet
    // Split the line into parts
    std::vector<std::string> elements = splitLine(line);
    std::string loadString = elements[loadIndex];
    if (loadString.find('.') != std::string::npos) {
      // Ensure the load value is valid
      try {
        double loadOfLine = std::stod(elements[loadIndex]);
        if (loadOfLine >= currentLoad) {
          // Found a line with a load >= current load, stop further processing
          foundLargerLoad = true;
          break;
        }
      } catch (const std::invalid_argument &e) {
        throw std::runtime_error("Invalid load value in line: " + line);
      }
    }
    // Else continue
  }

  inputFile.close(); // Close the input file

  if (!foundLargerLoad) {
    // If no match is found, there's nothing to trim, return early
    return;
  }

  // Reopen the file for writing, this will truncate the file
  std::ofstream outputFile(filePath, std::ios::trunc);
  if (!outputFile.is_open()) {
    throw std::runtime_error("Could not open file for writing");
  }

  // Write back the lines before the match to the file
  // C++ does not easily deal with deleting sections of a file, so for simpcity,
  // we re-write everything we want to keep.
  for (const auto &savedLine : lines) {
    outputFile << savedLine << "\n";
  }

  outputFile.close(); // Close the file after truncating
}

void saveConfigFile(Config conf) {
  // Construct the full file path
  fs::path filePath =
      fs::path(getOutputPath(conf.name, findOutputPath())) / CONFIGNAME;

  // Check if both paths exist before checking equivalence
  if (fs::exists(conf.configPath) && fs::exists(filePath)) {
    // Check if the source and destination paths are the same
    if (fs::equivalent(conf.configPath, filePath)) {
      return;
    }
  }

  try {
    std::ifstream src(conf.configPath, std::ios::binary);
    std::ofstream dst(filePath, std::ios::binary);

    if (!src) {
      if (conf.configPath.empty()) {
        std::cout << "No config path specified." << std::endl;
      } else {
        std::cerr << "Failed to open source file: " << conf.configPath
                  << std::endl;
      }
      return;
    }

    if (!dst) {
      std::cerr << "Failed to open destination file: " << filePath << std::endl;
      std::cerr << "Check if the output directory exists and you have "
                   "permission to write."
                << std::endl;
      return;
    }

    dst << src.rdbuf(); // Copy the content
  } catch (const std::exception &e) {
    std::cerr << "Exception occurred: " << e.what() << std::endl;
  }
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

// ChatGPT magic. Converts everything into strings
template <typename... Args>
std::vector<std::string> createStringVector(Args &&...args) {
  std::vector<std::string> vec;
  (vec.push_back([=] {
    std::ostringstream oss;
    // Set precision to 15 (you can adjust this as needed)
    oss.precision(13);
    oss << args;
    return oss.str();
  }()),
   ...);
  return vec;
}

// This writes any information we want to one line of the cvs file
void writeToCsv(std::ofstream &file, const Simulation &s) {
  auto lineData = getStringVector(s);
  writeLineToCsv(file, lineData);
}

std::vector<std::string> getStringVector(const Simulation &s) {
  // Must be in the same order as getCsvCols
  auto lineData = createStringVector(
      s.mesh.load, s.mesh.averageEnergy, s.mesh.delAvgEnergy, s.mesh.maxEnergy,
      s.mesh.averageRSS, s.mesh.nrPlasticChanges, s.FIRERep.nrIter,
      s.LBFGSRep.nrIter, s.CGRep.nrIter,                //
      s.FIRERep.nfev, s.LBFGSRep.nfev, s.CGRep.nfev,    //
      s.FIRERep.tType, s.LBFGSRep.tType, s.CGRep.tType, //
      s.timer.RTString(), s.timer.RTString("minimization", 7),
      s.timer.RTString("write", 7), s.timer.oldETRString, s.mesh.bounds[0],
      s.mesh.bounds[1], s.mesh.bounds[2], s.mesh.bounds[3], s.config.dtStart);
  return lineData;
}
std::vector<std::string> getCsvCols() {
  return {"Load",
          "Avg energy",
          "Avg energy change",
          "Max energy",
          "Avg RSS",
          "Nr plastic deformations",
          "Nr FIRE iterations",
          "Nr LBFGS iterations",
          "Nr CG iterations",
          "Nr FIRE func evals",
          "Nr LBFGS func evals",
          "Nr CG iterations",
          "FIRE Term reason",
          "LBFGS Term reason",
          "CG Term reason",
          "Run time",
          "Minimization time",
          "Write time",
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
bool insertHeaderIfNeeded(const std::string &filename) {
  std::ifstream fileIn(filename);

  if (!fileIn.is_open()) {
    // Create the file with only the header if it does not exist
    std::ofstream fileOut(filename);
    if (!fileOut.is_open()) {
      throw std::runtime_error("Unable to create file.");
    }
    writeCsvCols(fileOut);
    fileOut.close();
    return true;
  }

  // Read the first line from the file
  std::string firstLine;
  std::getline(fileIn, firstLine);
  fileIn.close();
  if (firstLine.empty()) {
    std::ofstream fileOut(filename);
    writeCsvCols(fileOut);
    fileOut.close();
    return true;
  }
  // Check if the header needs to be updated
  // Here we do a super lazy check. We just assume that if the first line starts
  // with the same character as the first character of the first column, things
  // are as they should be.
  else if (firstLine[0] != getCsvCols()[0][0]) {
    std::cout << "Warning! Modifying header in csv file: " << filename
              << std::endl;
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
    return true;
  } else {
    return false;
  }
}

void createCollection(const std::string folderPath,
                      const std::string destination,
                      const std::string extension,
                      const std::vector<double> &timestep) {

  std::vector<std::pair<int, fs::path>> filesWithNumbers;

  std::regex regexPattern(".*\\.([0-9]+)\\.vtu");

  for (const auto &entry : fs::directory_iterator(folderPath)) {
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

  fs::path relativePath =
      relative(folderPath,
               fs::path(folderPath)
                   .parent_path()); // Get the relative path from the parent

  std::ofstream outFile(destination + "/" + COLLECTIONNAME + ".pvd");
  outFile << "<?xml version=\"1.0\"?>\n";
  outFile << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
  outFile << "<Collection>\n";

  for (size_t i = 0; i < filesWithNumbers.size(); ++i) {
    double ts = timestep.size() > i ? timestep[i] : static_cast<double>(i);
    outFile << "<DataSet timestep=\"" << ts
            << "\" group=\"\" part=\"0\" file=\"" << DATAFOLDERPATH
            << filesWithNumbers[i].second.filename().string() << "\"/>\n";
  }

  outFile << "</Collection>\n";
  outFile << "</VTKFile>\n";
  outFile.close();
}

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