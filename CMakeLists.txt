cmake_minimum_required(VERSION 3.14)
project(MTS2D)

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
add_compile_options(-fopenmp)  # Enable OpenMP for parallel programming
add_compile_options(-fPIC)     # Position Independent Code for shared libraries

# Options specific to the Debug build
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0)      # No optimization, speeds up debugging
    add_compile_options(-g)       # Generate debug information
    add_compile_options(-Wall)    # Enable all warnings
    add_compile_options(-pthread) # Enable POSIX threads library
    
    # Additional tools for memory and behavior analysis on Linux
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        add_compile_options(-fsanitize=address)    # Memory access analyzer
        add_compile_options(-fsanitize=undefined)  # Undefined behavior analyzer
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address -fsanitize=undefined")
    endif()

# Options specific to the Release build
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-Ofast)           # Highest level of optimization, may break strict standards compliance
    add_compile_options(-march=native)    # Optimize code for the current machine's architecture
    add_compile_options(-funroll-loops)   # Unroll loops for speed
    add_compile_options(-floop-block)     # Reorder loops to improve cache utilization
    add_compile_options(-floop-interchange) # Interchange nested loops to improve cache usage
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -flto") # Enable Link Time Optimization
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto")
endif()

# Find packages
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "OpenMP not found")
endif()

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/alglib/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/FIREpp/FIRE
)
add_subdirectory(libs/cereal)
add_subdirectory(libs/eigen)
add_subdirectory(libs/itensorCmake)

# Source and header files
file(GLOB_RECURSE PROJECT_HEADERS src/**.h)
file(GLOB_RECURSE PROJECT_SOURCES src/**.cpp)

# Source grouping for IDEs
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

# stdc++fs is required for std::filesystem in C++17
# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE alglib itensor OpenMP::OpenMP_CXX stdc++fs cereal::cereal Eigen3::Eigen)
target_link_libraries(shared_objs PRIVATE alglib itensor OpenMP::OpenMP_CXX stdc++fs cereal::cereal Eigen3::Eigen)

# Testing
file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests/run)
target_link_libraries(test_${PROJECT_NAME} PRIVATE alglib itensor OpenMP::OpenMP_CXX stdc++fs cereal::cereal Eigen3::Eigen)