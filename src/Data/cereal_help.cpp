#include <cstdio> // For std::remove
#include <fstream>
#include <iostream>
#include <string>
#include <zlib.h>

constexpr size_t BUFFER_SIZE = 8192;

// **************************************************
// Save a file or raw string as a compressed .gz file
// **************************************************
bool saveCompressedGz(const std::string &filePath, const std::string &data,
                      bool isFile) {
  gzFile gzFile = gzopen(filePath.c_str(), "wb");
  if (!gzFile) {
    std::cerr << "Error: Could not open " << filePath << " for writing.\n";
    return false;
  }

  if (isFile) {
    std::ifstream inFile(data, std::ios::binary);
    if (!inFile) {
      std::cerr << "Error: Could not open " << data << " for reading.\n";
      gzclose(gzFile);
      return false;
    }

    char buffer[BUFFER_SIZE];
    while (inFile.read(buffer, BUFFER_SIZE) || inFile.gcount() > 0) {
      int written =
          gzwrite(gzFile, buffer, static_cast<unsigned>(inFile.gcount()));
      if (written <= 0) {
        int errnum;
        const char *errmsg = gzerror(gzFile, &errnum);
        std::cerr << "Error: gzwrite failed for " << filePath << ": " << errmsg
                  << "\n";
        gzclose(gzFile);
        return false;
      }
    }
    inFile.close();
  } else {
    if (gzwrite(gzFile, data.c_str(), static_cast<unsigned>(data.size())) ==
        0) {
      std::cerr << "Error: Failed to write compressed data to " << filePath
                << ".\n";
      gzclose(gzFile);
      return false;
    }
  }

  gzclose(gzFile);
  return true;
}

// **************************************************
// Load and decompress a .gz file into a string
// **************************************************
std::string loadCompressedGz(const std::string &filePath) {
  gzFile file = gzopen(filePath.c_str(), "rb");
  if (!file) {
    std::cerr << "Error: Could not open " << filePath << " for reading.\n";
    return "";
  }

  std::string data;
  char buffer[BUFFER_SIZE];
  int bytesRead;

  while ((bytesRead = gzread(file, buffer, BUFFER_SIZE)) > 0) {
    data.append(buffer, bytesRead);
  }

  if (bytesRead < 0) {
    std::cerr << "Error: Failed to read compressed file " << filePath << ".\n";
  }

  gzclose(file);
  return data;
}