#include "cereal_help.h"
#include <iostream>
#include <zlib.h>

// **************************************************
// Save XML data as a compressed .gz file
// **************************************************
void saveCompressedGz(const std::string &filePath, const std::string &xmlData) {
  // Open the output file using zlib's gzip functions
  gzFile gzFile = gzopen(filePath.c_str(), "wb");
  if (!gzFile) {
    std::cerr << "Error: Could not open " << filePath << " for writing.\n";
    return;
  }

  // Write the XML data to the compressed file
  if (gzwrite(gzFile, xmlData.c_str(), xmlData.size()) == 0) {
    std::cerr << "Error: Failed to write compressed data to " << filePath
              << "\n";
  }

  // Close the compressed file
  gzclose(gzFile);
}

// **************************************************
// Load and decompress a .gz XML file into a string
// **************************************************
std::string loadCompressedGz(const std::string &filePath) {
  // Open the .gz file for reading
  gzFile file = gzopen(filePath.c_str(), "rb");
  if (!file) {
    std::cerr << "Error: Could not open " << filePath << " for reading.\n";
    return "";
  }

  // Read decompressed data into a string
  std::string xmlData;
  char buffer[4096];
  int bytesRead;

  while ((bytesRead = gzread(file, buffer, sizeof(buffer))) > 0) {
    xmlData.append(buffer, bytesRead);
  }

  // Close the file
  gzclose(file);

  return xmlData;
}