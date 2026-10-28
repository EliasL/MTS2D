# Isolate third-party numerical code from the simulator's fast-math/LTO flags.
set_directory_properties(PROPERTIES COMPILE_OPTIONS "")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS "")
set(CMAKE_CXX_FLAGS "")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(BUILD_SHARED_LIBS OFF)
set(BUILD_STATIC_LIBS ON)
set(SUITESPARSE_ENABLE_PROJECTS "cholmod")
set(SUITESPARSE_USE_CUDA OFF)
set(SUITESPARSE_USE_FORTRAN OFF)
set(SUITESPARSE_USE_PYTHON OFF)
set(SUITESPARSE_USE_OPENMP OFF) # Parallel BLAS remains available.
set(SUITESPARSE_DEMOS OFF)

# Upstream configures headers in its source tree. Give each build its own copy
# so configurations cannot modify vendored sources or each other's headers.
set(suitesparse_source "${CMAKE_CURRENT_BINARY_DIR}/source")
file(COPY "${PROJECT_SOURCE_DIR}/libs/suitesparse/" DESTINATION "${suitesparse_source}")
add_subdirectory("${suitesparse_source}" "${CMAKE_CURRENT_BINARY_DIR}/build" EXCLUDE_FROM_ALL)
