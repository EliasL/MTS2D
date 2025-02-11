// ********************************************************************
// This header defines helper macros and a function to reduce duplicate
// writing of field names in serialization and deserialization.
// ********************************************************************
#ifndef CEREAL_HELPERS_H
#define CEREAL_HELPERS_H

#include <cereal/archives/xml.hpp>
#include <cereal/cereal.hpp>
#include <fstream> // Required for std::ifstream and std::ofstream
#include <iostream>
#include <sstream>
#include <string>

// This macro creates a name-value pair by converting the field name to a
// string. It avoids writing the field name twice (as both key and variable).
#define MAKE_NVP(field) cereal::make_nvp(#field, field)

// This macro calls loadWithDefault by automatically using the variable name
// as the key. It takes the archive, the field name, and its default value.
// It only does this if the cerialization function is used for loading though,
// not for saving.
#define LOAD_WITH_DEFAULT(ar, field, defaultValue)                             \
  do {                                                                         \
    if constexpr (cereal::traits::is_input_serializable<                       \
                      decltype(field), decltype(ar)>::value) {                 \
      loadWithDefault(ar, #field, field, defaultValue);                        \
    } else {                                                                   \
      ar(MAKE_NVP(field));                                                     \
    }                                                                          \
  } while (0)

// ***************************************
// Utility: Load field with default value
// ***************************************
template <class Archive, typename T>
void loadWithDefault(Archive &ar, const char *name, T &value,
                     const T &defaultValue) {
  T tempValue = defaultValue;
  try {
    ar(cereal::make_nvp(name, tempValue));
  } catch (const cereal::Exception &) {
    // Field missing: use default value
  }
  value = tempValue;
}

// ***************************************
// Explicit template instantiation (Fixes linker errors)
// ***************************************
// template std::string serializeToXml<Simulation>(const Simulation &);
// template void deserializeFromXml<Simulation>(const std::string &, Simulation
// &);

// ***************************************
// Serialize an object to XML
// ***************************************
template <typename T> std::string serializeToXml(const T &obj) {
  std::ostringstream xmlStream;
  {
    cereal::XMLOutputArchive oarchive(xmlStream);
    oarchive(cereal::make_nvp("Simulation", obj));
  }
  return xmlStream.str();
}

// ***************************************
// Deserialize an XML string into an object
// ***************************************
template <typename T>
void deserializeFromXml(const std::string &xmlData, T &obj) {
  std::istringstream iss(xmlData);
  {
    cereal::XMLInputArchive iarchive(iss);
    iarchive(cereal::make_nvp("Simulation", obj));
  }
}

// Function to save a compressed XML file inside .tar.gz
void saveCompressedGz(const std::string &filePath, const std::string &xmlData);

// Function to extract XML from a .tar.gz archive
std::string loadCompressedGz(const std::string &gzFile);

// Function to save data as either .xml or .gz
template <typename T>
void saveToFile(const std::string &filePath, const T &obj) {
  std::string xmlData = serializeToXml(obj);

  if (filePath.size() >= 3 && filePath.substr(filePath.size() - 3) == ".gz") {
    saveCompressedGz(filePath, xmlData);
  } else {
    std::ofstream outFile(filePath);
    if (!outFile) {
      std::cerr << "Error: Could not open " << filePath << " for writing.\n";
      return;
    }
    outFile << xmlData;
  }
}

// Function to load data from .xml or .gz, automatically decompressing if needed
template <typename T> void loadFromFile(const std::string &filePath, T &obj) {
  std::string xmlData;

  if (filePath.size() >= 3 && filePath.substr(filePath.size() - 3) == ".gz") {
    xmlData = loadCompressedGz(filePath);
  } else {
    std::ifstream inFile(filePath);
    if (!inFile) {
      std::cerr << "Error: Could not open " << filePath << " for reading.\n";
      return;
    }
    xmlData.assign((std::istreambuf_iterator<char>(inFile)),
                   std::istreambuf_iterator<char>());
  }

  deserializeFromXml(xmlData, obj);
}

#endif // CEREAL_HELPERS_H