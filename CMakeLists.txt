cmake_minimum_required(VERSION 3.5)
project(CrystalSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# Prevent easylogging from creating a default log file
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DELPP_NO_DEFAULT_LOG_FILE")

# Compiler options
# use -Ofast for performance builds
# use -Og for debugging and development=                    
add_compile_options(-Og -ggdb -fPIC)
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/libs/easylogging++
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
add_library(easyloggingpp STATIC libs/easylogging++/easylogging++.cpp)

# Find OpenMP
find_package(OpenMP)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE alglib easyloggingpp OpenMP::OpenMP_CXX)


# Testing
file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)

# Add the src directory to the include paths for your tests executable
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/tests/run)
target_link_libraries(test_${PROJECT_NAME} PRIVATE alglib easyloggingpp)
