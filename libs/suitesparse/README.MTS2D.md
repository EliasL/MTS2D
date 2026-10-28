# Bundled CHOLMOD

Upstream: https://github.com/DrTimothyAldenDavis/SuiteSparse/tree/v7.14.1

Version: SuiteSparse 7.14.1, including CHOLMOD 5.3.5.

This directory contains the unmodified upstream root CMakeLists.txt and license,
plus CHOLMOD, SuiteSparse_config, AMD, COLAMD, CAMD and CCOLAMD. CHOLMOD includes
its own SuiteSparse_metis ordering implementation. Package licenses and notices
are retained in their original directories; see LICENSE.txt and CHOLMOD/Doc.
No BLAS or LAPACK implementation is bundled.

`cmake/cholmod/CMakeLists.txt` builds the static CPU libraries with upstream
CMake. It disables demos, CUDA, Fortran compilation and internal OpenMP; the
external BLAS can still use threads. Library sources are optimized separately
from the simulator, without its fast-math/LTO flags, even in Debug builds.

Upstream generates configuration headers in its source directories. Each build
therefore copies this subset into its own ignored build directory before
configuration. This keeps the vendored files unchanged and prevents independent
builds from overwriting each other's headers. No network access is needed.

The upstream `v7.14.1.tar.gz` archive SHA256 is
`81e560e1f74546df139edb765b3f5bc865866da23062312ffe8fd821063c8397`.
