#ifndef COMPARE_MACROS_H
#define COMPARE_MACROS_H

#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>

namespace compare_detail {
template <typename T> inline std::string formatValue(const T &value) {
  if constexpr (std::is_same_v<T, bool>) {
    return value ? "true" : "false";
  } else if constexpr (std::is_floating_point_v<T>) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << value;
    return oss.str();
  } else if constexpr (std::is_integral_v<T>) {
    return std::to_string(value);
  } else if constexpr (std::is_enum_v<T>) {
    using U = std::underlying_type_t<T>;
    return std::to_string(static_cast<U>(value));
  } else {
    return std::string();
  }
}

template <typename T>
inline void appendValueDiff(std::string *debugMsg, int tabNumber,
                            const char *fieldName, const T &lhsVal,
                            const T &rhsVal) {
  if (!debugMsg) {
    return;
  }
  if constexpr (std::is_floating_point_v<T> || std::is_integral_v<T> ||
                std::is_enum_v<T> || std::is_same_v<T, bool>) {
    *debugMsg += std::string(tabNumber, '\t') + fieldName + " differs (lhs=" +
                 formatValue(lhsVal) + ", rhs=" + formatValue(rhsVal) +
                 "); \n";
  } else {
    *debugMsg += std::string(tabNumber, '\t') + fieldName + " differs; \n";
  }
}
} // namespace compare_detail

// ********************************************************************
// Macro to compare a field from lhs and rhs.
// This macro sets the flag 'equal' to false when the fields differ,
// and appends a debug message if a non-null debugMsg pointer is provided.
#define COMPARE_FIELD(field)                                                   \
  do {                                                                         \
    if (!(lhs.field == rhs.field)) {                                           \
      equal = false;                                                           \
      compare_detail::appendValueDiff(debugMsg, tabNumber, #field, lhs.field,  \
                                      rhs.field);                              \
    }                                                                          \
  } while (0)

#endif // COMPARE_MACROS_H
