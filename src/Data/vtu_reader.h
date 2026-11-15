#pragma once

#include <array>
#include <filesystem>
#include <string>
#include <vector>

// Geometry of a single-piece ASCII VTK unstructured grid.
struct VtuGrid {
  std::vector<std::array<double, 3>> points;
  std::vector<int> connectivity;
  std::vector<int> offsets;
  std::vector<int> types;
};

VtuGrid readVtu(const std::filesystem::path &path,
                const std::string &expectedSha256 = {});
