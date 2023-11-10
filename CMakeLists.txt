cmake_minimum_required(VERSION 3.0)
project(CrystalSimulation)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler options
add_compile_options(-O3 -ggdb -Xpreprocessor -fopenmp -fPIC)
include_directories(SYSTEM libs/alglib-cpp/src)

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
target_include_directories(shared_objs PRIVATE ${PROJECT_SOURCE_DIR}/src)

# Compile Alglib as a static library
file(GLOB ALGLIB_SOURCES libs/alglib-cpp/src/*.cpp)
add_library(alglib STATIC ${ALGLIB_SOURCES})

# Target Executable
add_executable(${PROJECT_NAME} src/main.cpp $<TARGET_OBJECTS:shared_objs>)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(${PROJECT_NAME} PRIVATE alglib -L/opt/homebrew/opt/llvm/lib/c++ -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++)
target_link_libraries(${PROJECT_NAME} PRIVATE -L/opt/homebrew/opt/libomp/lib -lomp)

# Testing
file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/tests/run)
target_link_libraries(test_${PROJECT_NAME} PRIVATE alglib)
