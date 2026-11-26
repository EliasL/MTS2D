#pragma once

#include <string_view>

inline bool usesTriangularReference(std::string_view energyFunction) {
  return energyFunction == "contiTriangular" || energyFunction == "triangularCubic";
}
