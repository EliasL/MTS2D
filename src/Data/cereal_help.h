// ********************************************************************
// This header defines helper macros and a function to reduce duplicate
// writing of field names in serialization and deserialization.
// ********************************************************************
#ifndef CEREAL_HELPERS_H
#define CEREAL_HELPERS_H

#include <cereal/cereal.hpp>

// This macro creates a name-value pair by converting the field name to a
// string. It avoids writing the field name twice (as both key and variable).
#define MAKE_NVP(field) cereal::make_nvp(#field, field)

// This macro calls loadWithDefault by automatically using the variable name
// as the key. It takes the archive, the field name, and its default value.
#define LOAD_WITH_DEFAULT(ar, field, defaultValue)                             \
  loadWithDefault(ar, #field, field, defaultValue)

// ********************************************************************
// This function loads a field from an archive and uses a default value if
// the field is missing. It creates a temporary variable for loading.
// ********************************************************************
template <class Archive, typename T>
void loadWithDefault(Archive &ar, const char *name, T &value,
                     const T &defaultValue) {
  // Create a temporary variable to hold the loaded value.
  T tempValue = defaultValue;

  // Try to load the value from the archive using its name.
  try {
    ar(cereal::make_nvp(name, tempValue));
  }
  // If the field is not found, the exception is caught and the default is kept.
  catch (const cereal::Exception &) {
    // Exception caught; tempValue remains as defaultValue.
  }

  // Assign the loaded (or default) value to the actual field.
  value = tempValue;
}

#endif // CEREAL_HELPERS_H