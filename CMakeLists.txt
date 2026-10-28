cmake_minimum_required(VERSION 3.14)
project(MTS2D)

# Verbose build information
message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Show compilers and versions
message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C compiler version: ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CXX compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

# Show OpenMP detection results
message(STATUS "OpenMP_CXX_FOUND: ${OpenMP_CXX_FOUND}")
message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")

# Fast math significantly speeds up calculations.
# The effect on accuracy is not fully known, but it seems minimal.
# The objective of this program is to find local minima. So long as the program
# converges to a minimum, small numerical inaccuracies are acceptable.
#
# Important: project code should not rely on NaN/Inf sentinels,
# std::isfinite/std::isnan, Eigen::allFinite(), or std::numeric_limits-based
# control flow. -ffast-math can invalidate those assumptions and make failures
# difficult to diagnose.
option(ENABLE_FAST_MATH "Use aggressive fast-math flags in Release builds" ON)
message(STATUS "ENABLE_FAST_MATH: ${ENABLE_FAST_MATH}")
if(ENABLE_FAST_MATH)
    set(FAST_MATH_FLAGS "-ffast-math")
else()
    set(FAST_MATH_FLAGS "")
endif()

option(ENABLE_NATIVE_ARCH "Optimize Release builds for the build host CPU" ON)
message(STATUS "ENABLE_NATIVE_ARCH: ${ENABLE_NATIVE_ARCH}")
if(ENABLE_NATIVE_ARCH)
    set(NATIVE_ARCH_FLAG "-march=native")
else()
    set(NATIVE_ARCH_FLAG "")
endif()

# Keep selected hot functions out of line so profilers can show meaningful
# stacks. Disable if you want maximum inlining and peak performance.
option(ENABLE_NOINLINE "Disable MTS_NOINLINE for max performance" ON)
message(STATUS "ENABLE_NOINLINE: ${ENABLE_NOINLINE}")
if(NOT ENABLE_NOINLINE)
    add_compile_definitions(MTS_NOINLINE=)
endif()

set(CMAKE_CXX_STANDARD 17)  
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


# Set options to avoid recompilation
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# Compiler options for all builds
#add_compile_options(-fopenmp)  # Enable OpenMP for parallel programming
add_compile_options(-fPIC)     # Position Independent Code for shared libraries

# Options specific to the Debug build
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0)      # No optimization, speeds up debugging
    add_compile_options(-g)       # Generate debug information
    add_compile_options(-Wall)    # Enable all warnings
    add_compile_options(-pthread) # Enable POSIX threads library
    
    add_compile_options(-fno-omit-frame-pointer) # Keep frame pointers for better stack traces
    #add_compile_options(-fsanitize=address)    # Memory access analyzer
    #add_compile_options(-fsanitize=undefined)  # Undefined behavior analyzer
    #add_compile_options(-fsanitize=thread)  # thread analyzer (does not work with openMP)
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined")
    #set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined")

# Options specific to the Release build
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    #-O3 
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG ${NATIVE_ARCH_FLAG} -flto=auto ${FAST_MATH_FLAGS} -funroll-loops -fstrict-aliasing -w")
    set(CMAKE_C_FLAGS_RELEASE   "-O3 -DNDEBUG ${NATIVE_ARCH_FLAG} -flto=auto ${FAST_MATH_FLAGS} -funroll-loops -fstrict-aliasing -w")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto=auto")

    # LTO at link time
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} -flto=auto")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto=auto")
endif()

# Find packages
# First try the default OpenMP detection.
find_package(OpenMP)

# On some macOS + Homebrew GCC setups, CMake may fail to detect OpenMP
# automatically. In that case, provide hints for libgomp and try again.
if(NOT OpenMP_CXX_FOUND)
    message(STATUS "Searching for OpenMP via Homebrew GCC/libgomp...")
    if(EXISTS "/opt/homebrew/opt/libomp")
        # Make this a list: each entry is one compiler option
        set(OpenMP_C_FLAGS "-Xpreprocessor;-fopenmp;-I/opt/homebrew/opt/libomp/include")
        set(OpenMP_C_LIB_NAMES "omp")

        set(OpenMP_CXX_FLAGS "-Xpreprocessor;-fopenmp;-I/opt/homebrew/opt/libomp/include")
        set(OpenMP_CXX_LIB_NAMES "omp")

        set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
        set(OpenMP_CXX_FOUND TRUE)
        set(OpenMP_C_FOUND TRUE)
        message(STATUS "Found OpenMP via Homebrew: /opt/homebrew/opt/libomp")
    endif()
endif()

# Create an OpenMP imported target if CMake didn't provide one
if(TARGET OpenMP::OpenMP_CXX)
    message(STATUS "Using existing OpenMP::OpenMP_CXX target.")
elseif(OpenMP_CXX_FOUND AND OpenMP_omp_LIBRARY)
    add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
        INTERFACE_COMPILE_OPTIONS "${OpenMP_CXX_FLAGS}"
        INTERFACE_LINK_LIBRARIES "${OpenMP_omp_LIBRARY}"
    )
    message(STATUS "Created OpenMP::OpenMP_CXX interface target (Homebrew libomp).")
else()
    message(STATUS "OpenMP not found: OpenMP::OpenMP_CXX will not be linked.")
endif()

# Use Zlib for compression
find_package(ZLIB REQUIRED)

option(MTS_ENABLE_CHOLMOD "Build bundled CHOLMOD (requires external BLAS/LAPACK)" ON)
if(MTS_ENABLE_CHOLMOD)
    if(CMAKE_VERSION VERSION_LESS 3.23)
        message(FATAL_ERROR "Bundled CHOLMOD requires CMake >= 3.23; use MTS_ENABLE_CHOLMOD=OFF for Eigen only.")
    endif()
    add_subdirectory(cmake/cholmod)
    set(tangent_backend_default CHOLMOD)
else()
    set(tangent_backend_default EIGEN)
endif()
set(MTS_TANGENT_BACKEND "${tangent_backend_default}" CACHE STRING "Default tangent backend: CHOLMOD or EIGEN")
set_property(CACHE MTS_TANGENT_BACKEND PROPERTY STRINGS CHOLMOD EIGEN)
if(NOT MTS_TANGENT_BACKEND MATCHES "^(CHOLMOD|EIGEN)$")
    message(FATAL_ERROR "MTS_TANGENT_BACKEND must be CHOLMOD or EIGEN.")
endif()
if(MTS_TANGENT_BACKEND STREQUAL "CHOLMOD" AND NOT MTS_ENABLE_CHOLMOD)
    message(FATAL_ERROR "MTS_TANGENT_BACKEND=CHOLMOD requires MTS_ENABLE_CHOLMOD=ON.")
endif()
message(STATUS "Default tangent backend: ${MTS_TANGENT_BACKEND}")

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/alglib/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/FIREpp/FIRE
)
add_subdirectory(libs/cereal)
add_subdirectory(libs/eigen)
#add_subdirectory(libs/itensorCmake)

# Source and header files
file(GLOB_RECURSE PROJECT_HEADERS CONFIGURE_DEPENDS src/**.h)
file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS src/**.cpp)

# Source grouping for IDEs (Integrated development environments)
source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})

# Exclude main.cpp from shared objects
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*main\\.cpp$")
add_library(shared_objs OBJECT ${PROJECT_SOURCES})

# Robust geometry predicates require strict floating-point arithmetic. Keep this
# isolated from fast-math/LTO in the simulator's other translation units.
set_source_files_properties(src/Mesh/mesh_delaunay.cpp PROPERTIES
    COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-fno-lto"
)
set_property(SOURCE src/Mesh/mesh_delaunay.cpp APPEND PROPERTY
    INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/libs/cdt/include")

# Main executable
add_executable(${PROJECT_NAME} src/main.cpp $<TARGET_OBJECTS:shared_objs>)

# Compile Alglib as a static library
file(GLOB ALGLIB_SOURCES libs/alglib/src/*.cpp)
add_library(alglib STATIC ${ALGLIB_SOURCES})

# Disable warnings from library files
file(GLOB_RECURSE ALL_LIBRARY_SOURCES libs/*/*.cpp)
set_source_files_properties(${ALL_LIBRARY_SOURCES} PROPERTIES
    COMPILE_OPTIONS "-w;-Wno-class-memaccess"
)



# ***************************************
# Decide which filesystem library to link.
# On Clang (with libc++), std::filesystem is built in.
# On GCC, we need to link to stdc++fs.
# ***************************************
set(COMMON_LIBS 
    alglib
    OpenMP::OpenMP_CXX
    cereal::cereal
    Eigen3::Eigen
    ${ZLIB_LIBRARIES}
)
if(MTS_ENABLE_CHOLMOD)
    list(APPEND COMMON_LIBS SuiteSparse::CHOLMOD_static)
    target_compile_definitions(shared_objs PRIVATE MTS_HAVE_CHOLMOD=1)
endif()
if(MTS_TANGENT_BACKEND STREQUAL "CHOLMOD")
    target_compile_definitions(shared_objs PRIVATE MTS_DEFAULT_CHOLMOD=1)
endif()

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Not using Clang, linking stdc++fs.")
    list(APPEND COMMON_LIBS stdc++fs)
else()
    message(STATUS "Using Clang, so not linking stdc++fs.")
endif()

# Link libraries to targets
target_link_libraries(${PROJECT_NAME} PRIVATE ${COMMON_LIBS})
target_link_libraries(shared_objs PUBLIC ${COMMON_LIBS})
target_include_directories(shared_objs PUBLIC
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/libs/alglib/src
    ${PROJECT_SOURCE_DIR}/libs/FIREpp/FIRE
)

option(MTS_BUILD_TANGENT_BENCHMARK "Build the explicit tangent backend benchmark" OFF)
if(MTS_BUILD_TANGENT_BENCHMARK)
    add_executable(benchmark_tangent tools/benchmark_tangent.cpp $<TARGET_OBJECTS:shared_objs>)
    target_include_directories(benchmark_tangent PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_link_libraries(benchmark_tangent PRIVATE ${COMMON_LIBS})
endif()

# Testing — only build in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "tests/*.cpp" "tests/run/run_tests.cpp")
    add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
    target_link_libraries(test_${PROJECT_NAME} PRIVATE ${COMMON_LIBS})
    target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests/run)
    target_compile_definitions(test_${PROJECT_NAME} PRIVATE MTS2D_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
endif()
