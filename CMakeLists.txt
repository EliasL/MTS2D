cmake_minimum_required(VERSION 3.14)
project(MTS2D)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE "Debug")
endif()

#TODO https://www.alglib.net/translator/man/manual.cpp.html#gs_comm_smp_compiling
# and other optimizations

# Compiler options
add_compile_options(-fopenmp -fPIC)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-Ofast)
endif()

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/alglib/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/Eigen
)
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

# Include additional libraries
include(FetchContent)

# spdlog
FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.13.0 # specify the version you need
  )
  FetchContent_GetProperties(spdlog)
  if(NOT spdlog_POPULATED)
  FetchContent_Populate(spdlog)
  add_subdirectory(${spdlog_SOURCE_DIR} ${spdlog_BINARY_DIR})
  endif()


  # stdc++fs is required for std::filesystem
  # Link libraries
  target_link_libraries(${PROJECT_NAME} PRIVATE alglib OpenMP::OpenMP_CXX spdlog::spdlog stdc++fs)
  target_link_libraries(shared_objs PRIVATE alglib OpenMP::OpenMP_CXX spdlog::spdlog stdc++fs)
  
  # Testing
  file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
  add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
  target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests/run)
  target_link_libraries(test_${PROJECT_NAME} PRIVATE alglib OpenMP::OpenMP_CXX spdlog stdc++fs)
