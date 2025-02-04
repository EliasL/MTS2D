#ifndef COMPARE_MACROS_H
#define COMPARE_MACROS_H

// ********************************************************************
// Macro to compare a field from lhs and rhs.
// This macro sets the flag 'equal' to false when the fields differ,
// and appends a debug message if a non-null debugMsg pointer is provided.
#define COMPARE_FIELD(field)                                                   \
  do {                                                                         \
    if (!(lhs.field == rhs.field)) {                                           \
      equal = false;                                                           \
      if (debugMsg) {                                                          \
        *debugMsg += std::string(tabNumber, '\t') + #field " differs; \n";     \
      }                                                                        \
    }                                                                          \
  } while (0)

#endif // COMPARE_MACROS_H