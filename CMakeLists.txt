cmake_minimum_required(VERSION 3.11)
project(MTS2D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Debug")
endif()

# Compiler options
add_compile_options(-fopenmp -fPIC)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-Ofast -g)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-Ofast)
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/src SYSTEM libs/alglib/src)
# For memoization library
include_directories(${PROJECT_SOURCE_DIR}/libs/thread-safe-lru)

# Source and header files
file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cpp)

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

# Find OpenMP
find_package(OpenMP)

# Required for memoization
find_package(TBB)


# Include additional libraries
add_subdirectory(libs/indicators)
add_subdirectory(libs/spdlog)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE alglib OpenMP::OpenMP_CXX TBB spdlog::spdlog indicators)
target_link_libraries(shared_objs PRIVATE alglib OpenMP::OpenMP_CXX TBB spdlog::spdlog indicators)

# Testing
file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests/run)
target_link_libraries(test_${PROJECT_NAME} PRIVATE alglib OpenMP::OpenMP_CXX TBB spdlog)
