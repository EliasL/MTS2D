#ifndef DEBUG_COMPARE_H
#define DEBUG_COMPARE_H

#include "mesh.h"
#include "node.h"
#include "tElement.h"
#include <cstddef>
#include <string>

std::string debugCompare(const Node &lhs, const Node &rhs, int tabNumber = 0);
std::string debugCompare(const GhostNode &lhs, const GhostNode &rhs,
                         int tabNumber = 0);
std::string debugCompare(const TElement &lhs, const TElement &rhs,
                         int tabNumber = 0);
std::string debugCompare(const Mesh &lhs, const Mesh &rhs,
                         std::size_t maxLines = 100);

#endif
