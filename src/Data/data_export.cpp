#include "data_export.h"

#include "Data/param_parser.h"
#include "Mesh/mesh.h"
#include "settings.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ostream>
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
      "/data/elundheim/",  // PMMH
      "/data2/elundheim/", // PMMH
      "/Users/elias/Work/PhD/Code/localData/",
      "/lustre/fswork/projects/rech/bph/uog82gz/", // JeanZay
  };
  if (const char *home = std::getenv("HOME")) {
    paths.emplace_back(fs::path(home) / "simulation"); // Magi persistent home
  }

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
