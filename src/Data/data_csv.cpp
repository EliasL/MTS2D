#include "data_export.h"

#include "../Simulation/simulation.h"
#include "Data/param_parser.h"
#include "settings.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

std::string getFolderPath(const std::string &name, const std::string &dataPath,
                          const fs::path &subfolder);
void createBackupOfFile(const fs::path &file, const fs::path &backupDir,
                        std::size_t minSize);

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

// Helper function to split a line by commas
std::vector<std::string> splitLine(const std::string &line);
std::vector<std::string>
getStringVector(const Simulation &s, const std::vector<std::string> &headers);

static bool csvContainsLine(const std::string &filePath,
                            const std::string &target) {
  std::ifstream file(filePath);
  if (!file.is_open()) {
    return false;
  }
  std::string line;
  while (std::getline(file, line)) {
    if (line == target) {
      return true;
    }
  }
  return false;
}

static void appendHeaderChangeCommentIfNeeded(const std::string &filePath,
                                              const Simulation &s) {
  std::ifstream file(filePath);
  if (!file.is_open()) {
    return;
  }
  std::string header;
  if (!std::getline(file, header)) {
    return;
  }
  if (!header.empty() && header[0] == '#') {
    return;
  }
  const size_t fileHeaderCount = splitLine(header).size();
  const auto expectedHeaders = getCsvHeaders(s);
  if (fileHeaderCount == expectedHeaders.size()) {
    return;
  }

  const std::string comment =
      std::string("#HEADER:") + join(expectedHeaders, ",");
  if (csvContainsLine(filePath, comment)) {
    return;
  }

  std::ofstream out(filePath, std::ios::app);
  if (!out.is_open()) {
    return;
  }
  if (!isQuiet()) {
    std::cout << "CSV header size changed (" << fileHeaderCount << " -> "
              << expectedHeaders.size() << "); adding comment header."
              << std::endl;
  }
  out << comment << "\n";
}

static fs::path macroCsvPath(const std::string &folderName,
                             const std::string &dataPath,
                             const std::string &subFolder) {
  return fs::path(getFolderPath(folderName, dataPath, subFolder)) /
         (MACRODATANAME + std::string(".csv"));
}

// Function to initialize a CSV file for writing
std::ofstream initCsvFile(const std::string &folderName,
                          const std::string &dataPath, const Simulation &s,
                          const std::string subFolder,
                          bool appendHeaderChangeComment) {
  const fs::path filePath = macroCsvPath(folderName, dataPath, subFolder);
  bool headerWasWritten = insertHeaderIfNeeded(filePath, s);

  if (!headerWasWritten) {
    // If we start from a dump, we need to trim the csv file
    // Before we do, we will create a backup of the file
    fs::path backupDir = getBackupPath(s.simName, dataPath);
    // create backup if it is larger than 100KB
    createBackupOfFile(filePath, backupDir, 100 * 1024); // 100KB
    trimCsvFile(filePath, s);
    if (appendHeaderChangeComment) {
      appendHeaderChangeCommentIfNeeded(filePath, s);
    }
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

std::vector<std::string> readCsvHeaders(const std::string &folderName,
                                        const std::string &dataPath,
                                        const std::string subFolder) {
  const fs::path filePath = macroCsvPath(folderName, dataPath, subFolder);
  std::ifstream file(filePath);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open CSV file: " + filePath.string());
  }

  std::string line;
  while (std::getline(file, line)) {
    if (!line.empty() && line[0] == '#') {
      continue;
    }
    if (!line.empty()) {
      return splitLine(line);
    }
  }

  throw std::runtime_error("CSV file has no header: " + filePath.string());
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
  size_t loadStepIndex = std::string::npos;

  // Read line by line
  while (std::getline(inputFile, line)) {
    if (!line.empty() && line[0] == '#') {
      lines.push_back(line);
      continue;
    }
    // Split the line into columns
    std::vector<std::string> elements = splitLine(line);

    if (firstLine) {
      auto it = std::find(elements.begin(), elements.end(), "load_step");
      if (it == elements.end()) {
        throw std::runtime_error(
            "CSV header missing required load_step column: " + filePath);
      }
      loadStepIndex = static_cast<size_t>(std::distance(elements.begin(), it));
      lines.push_back(line);
      firstLine = false;
      continue;
    }

    try {
      // Convert the relevant column to an integer
      if (loadStepIndex >= elements.size()) {
        throw std::runtime_error("Missing load_step value in line: " + line);
      }
      const std::string &loadStepStr = elements[loadStepIndex];
      if (loadStepStr.empty()) {
        throw std::runtime_error("Missing load_step value in line: " + line);
      }
      int lineLoadStep = std::stoi(loadStepStr);

      // If the line's loadStep is larger than the current loadStep,
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
    lines.push_back(line);
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

// Function to write line to CSV using an open file stream
void writeLineToCsv(std::ofstream &file,
                    const std::vector<std::string> &strings) {
  if (!file.is_open()) {
    throw std::runtime_error("File stream is not open.");
  }

  std::string line = join(strings, ",");
  file << line << '\n';
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

void writeToCsv(std::ofstream &file, const Simulation &s) {
  const auto lineData = getStringVector(s);
  writeLineToCsv(file, lineData);
}

void writeToCsv(std::ofstream &file, const Simulation &s,
                const std::vector<std::string> &headers) {
  const auto lineData = getStringVector(s, headers);
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

std::vector<std::string>
getStringVector(const Simulation &s, const std::vector<std::string> &headers) {
  const auto &cols = s.getCsvColumns();
  std::vector<std::string> row;
  row.reserve(headers.size());
  for (const auto &header : headers) {
    auto it = std::find_if(cols.begin(), cols.end(), [&](const CsvColumn &col) {
      return col.name == header;
    });
    if (it == cols.end()) {
      throw std::runtime_error("CSV header column is not available: " + header);
    }
    row.push_back(it->getter(s));
  }
  return row;
}

bool insertHeaderIfNeeded(const std::string &filename, const Simulation &s) {
  const fs::path p = fs::path(filename);

  if (fs::exists(p)) {
    std::error_code ec;
    const auto sz = fs::file_size(p, ec);
    if (!ec && sz > 0) {
      return false;
    }
  }

  std::ofstream fileOut(filename, std::ios::trunc);
  if (!fileOut.is_open()) {
    throw std::runtime_error("Unable to create/open file with header: " +
                             filename);
  }

  writeCsvHeaders(fileOut, s);
  return true;
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
