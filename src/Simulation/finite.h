#pragma once
#include <cstddef>
#include <cstdint>

// std::isfinite can be optimized away under Release fast-math flags.
inline bool finiteValue(double value) {
  static_assert(sizeof(double) == sizeof(std::uint64_t));
  volatile double storedValue = value;
  std::uint64_t bits = 0;
  auto *destination = reinterpret_cast<unsigned char *>(&bits);
  const auto *bytes =
      reinterpret_cast<const volatile unsigned char *>(&storedValue);
  for (std::size_t i = 0; i < sizeof(value); ++i) {
    destination[i] = bytes[i];
  }
  constexpr std::uint64_t exponentMask = UINT64_C(0x7ff0000000000000);
  return (bits & exponentMask) != exponentMask;
}
