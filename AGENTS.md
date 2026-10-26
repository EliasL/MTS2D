# Project instructions

- Only add high-value tests that can catch plausible, consequential errors; do
  not add tests merely for coverage or trivial implementation details. Prefer
  regressions for real bugs. The simulator's scientific/numerical code warrants
  more careful testing: preserve its existing tests, unlike the lighter test
  policy used for management and plotting.

- Do not use `std::numeric_limits`: it is incompatible with this project's
  optimization level. Use the C floating-point limits from `<cfloat>` (for
  example, `DBL_MAX`) instead.
