# Project instructions

- MTS2D executes a prepared run folder and owns scientific output, checkpoint
  recovery and numerical stopping criteria. MTSManager owns job preparation,
  retained executables and scheduler records.
- Coordinate configuration, CLI and output changes with the manager schema and
  plotting readers. Keep saved configurations and checkpoint identity stable.
  Use the native Debug test preset for scientific checks.

- Only add high-value tests that can catch plausible, consequential errors; do
  not add tests merely for coverage or trivial implementation details. Prefer
  regressions for real bugs. The simulator's scientific/numerical code warrants
  more careful testing: preserve its existing tests, unlike the lighter test
  policy used for management and plotting.

- Do not use `std::numeric_limits`: it is incompatible with this project's
  optimization level. Use the C floating-point limits from `<cfloat>` (for
  example, `DBL_MAX`) instead.

# Element tracking

Follow [Element identity and transport](../AGENTS.md#element-identity-and-transport).
Preserve parent-to-child `eIndex` correspondence when changing edge-flip
reconnection or VTU export. Long/short transport changes the geometry assigned
to a history, not whether that history retains its element identity.
