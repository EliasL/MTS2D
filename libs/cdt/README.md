# Bundled CDT

Source: https://github.com/artem-ogre/CDT

Pinned commit: `15c15405329e50303481df8f33e1bb6764e30f74`.
Only the upstream `CDT/include/` files and their licenses are vendored,
without modifications. See `LICENSE`, `LICENSE-LLVM.txt`, and
`LICENSE-predicates.txt` for the respective terms.

MTS2D uses ordinary Delaunay triangulation (`insertVertices` followed by
`eraseSuperTriangle`). It does not use constraints, conforming triangulation,
or refinement, and does not add physical nodes. Periodic halo generation and
validation remain in `src/Mesh/mesh_delaunay.cpp`.

That translation unit is compiled without fast-math, implicit floating-point
contraction, or LTO. This protects the robust predicates and upstream numeric
constants from the simulator's aggressive Release optimization flags.
