# Tangent backend comparison

The optional `benchmark_tangent` target compares both backends in one Release
executable on the same assembled periodic tangent and right-hand side. It checks
physical solution agreement and the residual of the original linear system;
transformed coordinates are not compared because the permutations differ.

The measurements below were taken on macOS arm64 with Apple Clang 17 and the
system Accelerate BLAS/LAPACK. SuiteSparse 7.14.1 supplies CHOLMOD 5.3.5.
The simulator uses its normal Release flags; SuiteSparse uses `-O3 -DNDEBUG`
without fast-math/LTO. `VECLIB_MAXIMUM_THREADS`, `OPENBLAS_NUM_THREADS` and
`OMP_NUM_THREADS` were set to `1`. Each entry is the median of three runs;
each run times 20 pairs of coordinate/gradient half-solves.

| Mesh | Active DOFs | Eigen analysis + factor | CHOLMOD analysis + factor | Eigen half-solve pair | CHOLMOD half-solve pair |
|---|---:|---:|---:|---:|---:|
| 20 x 20 | 798 | 0.887 ms | 1.086 ms | 0.061 ms | 0.030 ms |
| 100 x 100 | 19,998 | 120.596 ms | 46.526 ms | 2.908 ms | 1.573 ms |
| 200 x 200 | 79,998 | 928.435 ms | 220.537 ms | 14.172 ms | 7.091 ms |

[Raw timings](tangent_backend_comparison.csv) separate assembly, analysis and
numerical factorization. At 200 x 200, numerical factorization alone was
911.815 ms for Eigen and 44.466 ms for CHOLMOD. CHOLMOD's more expensive
AMD/METIS ordering analysis reduces the benefit for a fresh preparation to
about 4.2 times. Eigen is slightly faster to prepare the smallest matrix.
Relative residuals ranged from approximately 1.3e-15 to 3.6e-15.
These are tangent-operation measurements, not whole-simulation speedups.

Both the CHOLMOD-default and Eigen-only builds passed all 116 scientific tests,
including coordinate round trips, physical inverse agreement, transformed
energy derivatives and diagonal-shift recovery. The added backend comparisons
use small matrices and do not run performance benchmarks in the test suite.

A separate Release simulation comparison used a 20 x 20 periodic mesh, seed 0,
QDSD 0.01, initial guess noise 0.001, no reconnection, and five shear increments
of 0.02. Both backends reported identical total energies and average stress
components at CSV precision in all six states. Both met the initialization
force tolerance of 1e-6 and subsequent tolerance of 1e-7, without preconditioner
fallbacks. Different results near plastic instabilities remain possible from
floating-point differences, so this elastic comparison does not establish
identical long-run avalanche trajectories.
