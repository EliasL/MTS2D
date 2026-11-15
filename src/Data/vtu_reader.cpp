#include "vtu_reader.h"
#include <cereal/macros.hpp>
#include <cereal/external/rapidxml/rapidxml.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <locale>
#include <openssl/evp.h>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace {
namespace xml = cereal::rapidxml;
using Node = xml::xml_node<>;

std::string_view attribute(const Node &node, const char *name) {
  const auto *value = node.first_attribute(name);
  return value ? value->value() : "";
}

const Node &child(const Node &node, const char *name) {
  const auto *result = node.first_node(name);
  if (!result || result->next_sibling(name)) {
    throw std::invalid_argument(std::string("VTU requires one ") + name + " section.");
  }
  return *result;
}

const Node &cellArray(const Node &cells, std::string_view name) {
  const Node *result = nullptr;
  for (auto *array = cells.first_node("DataArray"); array; array = array->next_sibling("DataArray")) {
    if (attribute(*array, "Name") != name) continue;
    if (result) throw std::invalid_argument("Duplicate VTU cell array: " + std::string(name));
    result = array;
  }
  if (!result) throw std::invalid_argument("Missing VTU cell array: " + std::string(name));
  return *result;
}

template <typename T> std::vector<T> values(const Node &array) {
  if (attribute(array, "format") != "ascii") {
    throw std::invalid_argument("VTU geometry must use ASCII DataArrays.");
  }
  std::istringstream input(array.value());
  input.imbue(std::locale::classic());
  std::vector<T> result;
  while (input >> std::ws && !input.eof()) {
    T value;
    if (!(input >> value) || !std::isfinite(value)) {
      throw std::invalid_argument("Invalid number in VTU geometry.");
    }
    result.push_back(value);
  }
  return result;
}

size_t count(const Node &piece, const char *name) {
  std::istringstream input{std::string(attribute(piece, name))};
  int number;
  if (!(input >> number) || number < 0 || !(input >> std::ws).eof()) {
    throw std::invalid_argument(std::string("Invalid VTU ") + name + ".");
  }
  return static_cast<size_t>(number);
}

std::string sha256(const std::string &bytes) {
  std::array<unsigned char, EVP_MAX_MD_SIZE> digest{};
  unsigned int size = 0;
  if (!EVP_Digest(bytes.data(), bytes.size(), digest.data(), &size, EVP_sha256(), nullptr)) {
    throw std::runtime_error("Could not calculate VTU checksum.");
  }
  std::ostringstream hex;
  hex << std::hex << std::setfill('0');
  for (unsigned int i = 0; i < size; ++i) hex << std::setw(2) << static_cast<int>(digest[i]);
  return hex.str();
}
} // namespace

VtuGrid readVtu(const std::filesystem::path &path, const std::string &expectedSha256) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::invalid_argument("Cannot open VTU file: " + path.string());
  std::string bytes{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
  if (!expectedSha256.empty() && sha256(bytes) != expectedSha256) {
    throw std::invalid_argument("The initialMesh.vtu checksum differs from config.conf.");
  }
  if (bytes.find('\0') != std::string::npos) throw std::invalid_argument("VTU must be ASCII XML.");
  xml::xml_document<> document;
  try {
    document.parse<xml::parse_validate_closing_tags>(bytes.data());
  } catch (const xml::parse_error &error) {
    throw std::invalid_argument(std::string("Invalid VTU XML: ") + error.what());
  }
  const auto &file = child(document, "VTKFile");
  if (attribute(file, "type") != "UnstructuredGrid") {
    throw std::invalid_argument("VTU must contain an UnstructuredGrid.");
  }
  const auto &piece = child(child(file, "UnstructuredGrid"), "Piece");
  const auto &points = child(child(piece, "Points"), "DataArray");
  const auto &cells = child(piece, "Cells");
  const auto coordinates = values<double>(points);
  const size_t pointCount = count(piece, "NumberOfPoints");
  const size_t cellCount = count(piece, "NumberOfCells");
  if (attribute(points, "NumberOfComponents") != "3" || coordinates.size() != 3 * pointCount) {
    throw std::invalid_argument("VTU requires three coordinates per point.");
  }
  VtuGrid grid;
  grid.points.reserve(pointCount);
  for (size_t i = 0; i < coordinates.size(); i += 3) {
    grid.points.push_back({coordinates[i], coordinates[i + 1], coordinates[i + 2]});
  }
  grid.connectivity = values<int>(cellArray(cells, "connectivity"));
  grid.offsets = values<int>(cellArray(cells, "offsets"));
  grid.types = values<int>(cellArray(cells, "types"));
  if (grid.offsets.size() != cellCount || grid.types.size() != cellCount) {
    throw std::invalid_argument("VTU cell counts do not match the cell arrays.");
  }
  int previous = 0;
  for (int offset : grid.offsets) {
    if (offset <= previous || static_cast<size_t>(offset) > grid.connectivity.size()) {
      throw std::invalid_argument("Invalid VTU cell offset.");
    }
    previous = offset;
  }
  if (static_cast<size_t>(previous) != grid.connectivity.size()) {
    throw std::invalid_argument("VTU offsets do not cover the connectivity array.");
  }
  for (int index : grid.connectivity) {
    if (index < 0 || static_cast<size_t>(index) >= pointCount) {
      throw std::invalid_argument("VTU connectivity references a missing point.");
    }
  }
  return grid;
}
