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

# Target Executable
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE -L/opt/homebrew/opt/llvm/lib/c++ -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++)
target_link_libraries(${PROJECT_NAME} PRIVATE -L/opt/homebrew/opt/libomp/lib -lomp)

# Link libraries
file(GLOB ALGLIB_SOURCES libs/alglib-cpp/src/*.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${ALGLIB_SOURCES})
