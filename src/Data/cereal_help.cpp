#include "cereal_help.h"
#include <iostream>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <string>
#include <zlib.h>

std::string compressFile(const std::string &filePath) {

  std::ifstream inFile(filePath, std::ios::binary);
  if (!inFile) {
    std::cerr << "Error: Could not open " << filePath << " for reading.\n";
    return "";
  }

  std::string outFilePath = filePath + ".gz";
  gzFile gzOut = gzopen(outFilePath.c_str(), "wb");
  if (!gzOut) {
    std::cerr << "Error: Could not open " << outFilePath << " for writing.\n";
    return "";
  }

  constexpr size_t BUFFER_SIZE = 8192;
  char buffer[BUFFER_SIZE];

  while (inFile.read(buffer, BUFFER_SIZE) || inFile.gcount() > 0) {
    if (gzwrite(gzOut, buffer, inFile.gcount()) == 0) {
      std::cerr << "Error: Failed to write compressed data to " << outFilePath
                << "\n";
      gzclose(gzOut);
      return "";
    }
  }

  inFile.close();
  gzclose(gzOut);

  // Remove original file
  if (std::remove(filePath.c_str()) != 0) {
    std::cerr << "Warning: Could not remove original file " << filePath << "\n";
  }
  return outFilePath;
}

// **************************************************
// Save data as a compressed .gz file
// **************************************************
void saveCompressedGz(const std::string &filePath, const std::string &data) {
  // Open the output file using zlib's gzip functions
  gzFile gzFile = gzopen(filePath.c_str(), "wb");
  if (!gzFile) {
    std::cerr << "Error: Could not open " << filePath << " for writing.\n";
    return;
  }

  // Write the data to the compressed file
  if (gzwrite(gzFile, data.c_str(), data.size()) == 0) {
    std::cerr << "Error: Failed to write compressed data to " << filePath
              << "\n";
  }

  // Close the compressed file
  gzclose(gzFile);
}

// **************************************************
// Load and decompress a .gz file into a string
// **************************************************
std::string loadCompressedGz(const std::string &filePath) {
  // Open the .gz file for reading
  gzFile file = gzopen(filePath.c_str(), "rb");
  if (!file) {
    std::cerr << "Error: Could not open " << filePath << " for reading.\n";
    return "";
  }

  // Read decompressed data into a string
  std::string data;
  char buffer[4096];
  int bytesRead;

  while ((bytesRead = gzread(file, buffer, sizeof(buffer))) > 0) {
    data.append(buffer, bytesRead);
  }

  // Close the file
  gzclose(file);

  return data;
}