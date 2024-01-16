cmake_minimum_required(VERSION 3.11)
project(CrystalSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-Og -ggdb -fPIC -pg)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-Ofast -pg)
endif()


include_directories(
    ${PROJECT_SOURCE_DIR}/src
    SYSTEM libs/alglib-cpp/src
)
# Include directories
file(GLOB_RECURSE PROJECT_HEADERS src/*.h)
file(GLOB_RECURSE PROJECT_SOURCES src/*.cpp)
file(GLOB PROJECT_CONFIGS CMakeLists.txt
                         .gitattributes .gitignore .gitmodules)

source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})


# Create an OBJECT library for shared sources, excluding main.cpp
list(FILTER PROJECT_SOURCES EXCLUDE REGEX ".*main\\.cpp$")
add_library(shared_objs OBJECT ${PROJECT_SOURCES} ${PROJECT_HEADERS})


# Target Executable
add_executable(${PROJECT_NAME} src/main.cpp $<TARGET_OBJECTS:shared_objs>)

# Compile Alglib as a static library
file(GLOB ALGLIB_SOURCES libs/alglib-cpp/src/*.cpp)
add_library(alglib STATIC ${ALGLIB_SOURCES})

# Compile Easylogging++ as a static library
#add_library(easyloggingpp STATIC libs/easylogging++/easylogging++.cpp)

# Find OpenMP
find_package(OpenMP)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE alglib OpenMP::OpenMP_CXX)

include(FetchContent)

FetchContent_Declare(
  spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.13.0
)
FetchContent_Declare(
  indicators
  GIT_REPOSITORY https://github.com/p-ranav/indicators.git
  GIT_TAG v2.1 
)
FetchContent_MakeAvailable(spdlog indicators)

if(NOT TARGET spdlog)
  find_package(spdlog REQUIRED)
endif()
if(NOT TARGET indicators)
  find_package(indicators REQUIRED)
endif()

# Link spdlog to shared_objs
target_link_libraries(shared_objs PRIVATE spdlog::spdlog indicators)

# Link spdlog to the main executable
target_link_libraries(${PROJECT_NAME} PRIVATE spdlog::spdlog indicators)


# Testing
file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)

# Add the src directory to the include paths for your tests executable
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/tests/run)
target_link_libraries(test_${PROJECT_NAME} PRIVATE alglib spdlog)
