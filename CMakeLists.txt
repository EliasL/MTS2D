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
option(ENABLE_FAST_MATH "Use aggressive fast-math flags in Release builds" ON)
message(STATUS "ENABLE_FAST_MATH: ${ENABLE_FAST_MATH}")

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
    # Baseline: high optimisation, standards-compliant FP, NDEBUG, arch tuning, LTO
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -march=native -flto=auto")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -DNDEBUG -march=native -flto=auto")

    # Optional fast-math mode
    if(ENABLE_FAST_MATH)
        # Slightly less extreme than -Ofast but still very aggressive
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -ffast-math")
        # If you *really* want -Ofast instead, you could do:
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast")
        # set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Ofast")
    endif()

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

# Set CGAL path manually (since we're not using system installation)
set(CGAL_DIR ${PROJECT_SOURCE_DIR}/libs/cgal)
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
list(APPEND CMAKE_PREFIX_PATH ${CGAL_DIR})

# Include CGAL
find_package(CGAL REQUIRED)
    
# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/alglib/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/FIREpp/FIRE
    SYSTEM ${CGAL_DIR}/include
)
add_subdirectory(libs/cereal)
add_subdirectory(libs/eigen)
#add_subdirectory(libs/itensorCmake)

# Source and header files
file(GLOB_RECURSE PROJECT_HEADERS src/**.h)
file(GLOB_RECURSE PROJECT_SOURCES src/**.cpp)

# Source grouping for IDEs (Integrated development environments)
source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})

# Exclude main.cpp from shared objects
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*main\\.cpp$")
add_library(shared_objs OBJECT ${PROJECT_SOURCES})

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
    cereal::cereal
    Eigen3::Eigen
    Threads::Threads
    CGAL::CGAL
    ${ZLIB_LIBRARIES}
)

# Add OpenMP library depending on the toolchain.
if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    # We configured APPLE_OPENMP_LIB above to point at libomp.dylib.
    list(APPEND COMMON_LIBS ${APPLE_OPENMP_LIB})
else()
    list(APPEND COMMON_LIBS OpenMP::OpenMP_CXX)
endif()

if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "Not using Clang, linking stdc++fs.")
    list(APPEND COMMON_LIBS stdc++fs)
else()
    message(STATUS "Using Clang, so not linking stdc++fs.")
endif()

# Link libraries to targets
target_link_libraries(${PROJECT_NAME} PRIVATE ${COMMON_LIBS})
target_link_libraries(shared_objs PRIVATE ${COMMON_LIBS})

# Testing — only build in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
    add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
    target_link_libraries(test_${PROJECT_NAME} PRIVATE ${COMMON_LIBS})
    target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests/run)
endif()