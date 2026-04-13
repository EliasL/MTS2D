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
if(ENABLE_FAST_MATH)
    set(FAST_MATH_FLAGS "-ffast-math")
else()
    set(FAST_MATH_FLAGS "")
endif()

# Profile-guided optimization (PGO)
set(PGO_MODE "OFF" CACHE STRING "PGO mode: OFF, GEN, USE")
set_property(CACHE PGO_MODE PROPERTY STRINGS OFF GEN USE)
set(PGO_DIR "${CMAKE_BINARY_DIR}/pgo" CACHE PATH "PGO profile data directory")
message(STATUS "PGO_MODE: ${PGO_MODE}")
message(STATUS "PGO_DIR: ${PGO_DIR}")

string(TOUPPER "${PGO_MODE}" PGO_MODE_UPPER)
if(NOT PGO_MODE_UPPER STREQUAL "OFF" AND
   NOT PGO_MODE_UPPER STREQUAL "GEN" AND
   NOT PGO_MODE_UPPER STREQUAL "USE")
    message(FATAL_ERROR "Invalid PGO_MODE: ${PGO_MODE}. Use OFF, GEN, or USE.")
endif()

set(PGO_C_FLAGS "")
set(PGO_CXX_FLAGS "")
set(PGO_LINK_FLAGS "")

if(PGO_MODE_UPPER STREQUAL "GEN" OR PGO_MODE_UPPER STREQUAL "USE")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(PGO_RAW_FILE "${PGO_DIR}/default.profraw")
        set(PGO_PROFDATA_FILE "${PGO_DIR}/default.profdata")
        if(PGO_MODE_UPPER STREQUAL "GEN")
            set(PGO_C_FLAGS "-fprofile-instr-generate=${PGO_RAW_FILE}")
            set(PGO_CXX_FLAGS "-fprofile-instr-generate=${PGO_RAW_FILE}")
            set(PGO_LINK_FLAGS "-fprofile-instr-generate=${PGO_RAW_FILE}")
        else()
            set(PGO_C_FLAGS "-fprofile-instr-use=${PGO_PROFDATA_FILE}")
            set(PGO_CXX_FLAGS "-fprofile-instr-use=${PGO_PROFDATA_FILE}")
            set(PGO_LINK_FLAGS "-fprofile-instr-use=${PGO_PROFDATA_FILE}")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        if(PGO_MODE_UPPER STREQUAL "GEN")
            set(PGO_C_FLAGS "-fprofile-generate=${PGO_DIR}")
            set(PGO_CXX_FLAGS "-fprofile-generate=${PGO_DIR}")
            set(PGO_LINK_FLAGS "-fprofile-generate=${PGO_DIR}")
        else()
            set(PGO_C_FLAGS "-fprofile-use=${PGO_DIR} -fprofile-correction")
            set(PGO_CXX_FLAGS "-fprofile-use=${PGO_DIR} -fprofile-correction")
            set(PGO_LINK_FLAGS "-fprofile-use=${PGO_DIR} -fprofile-correction")
        endif()
    else()
        message(WARNING "PGO requested, but compiler ${CMAKE_CXX_COMPILER_ID} is not supported.")
    endif()
endif()

# Keep selected hot functions out of line so profilers can show meaningful
# stacks. Disable if you want maximum inlining and peak performance.
option(ENABLE_NOINLINE "Disable MTS_NOINLINE for max performance" ON)
message(STATUS "ENABLE_NOINLINE: ${ENABLE_NOINLINE}")
if(NOT ENABLE_NOINLINE)
    add_compile_definitions(MTS_NOINLINE=)
endif()

option(IDE_LIGHTWEIGHT "Skip CGAL to speed up IDE parsing" OFF)
message(STATUS "IDE_LIGHTWEIGHT: ${IDE_LIGHTWEIGHT}")
if(IDE_LIGHTWEIGHT)
    add_compile_definitions(IDE_LIGHTWEIGHT)
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
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -flto=auto ${FAST_MATH_FLAGS} -funroll-loops -fstrict-aliasing -w ${PGO_CXX_FLAGS}")
    set(CMAKE_C_FLAGS_RELEASE   "-O3 -DNDEBUG -march=native -flto=auto ${FAST_MATH_FLAGS} -funroll-loops -fstrict-aliasing -w ${PGO_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto=auto ${PGO_LINK_FLAGS}")

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
list(APPEND CMAKE_PREFIX_PATH ${CGAL_DIR})

# Include CGAL (skip in IDE_LIGHTWEIGHT)
if(NOT IDE_LIGHTWEIGHT)
    find_package(CGAL REQUIRED)
endif()
    
# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/alglib/src
    SYSTEM ${PROJECT_SOURCE_DIR}/libs/FIREpp/FIRE
)
if(NOT IDE_LIGHTWEIGHT)
    include_directories(SYSTEM ${CGAL_DIR}/include)
endif()
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
    OpenMP::OpenMP_CXX
    cereal::cereal
    Eigen3::Eigen
    ${ZLIB_LIBRARIES}
)
if(NOT IDE_LIGHTWEIGHT)
    list(APPEND COMMON_LIBS CGAL::CGAL)
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

# Optional PGO training run target (uses configs/pgo_50x50.conf).
set(PGO_CONFIG_FILE "${PROJECT_SOURCE_DIR}/configs/pgo_50x50.conf")
add_custom_target(pgo-run
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/pgo_run
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR}
            $<TARGET_FILE:${PROJECT_NAME}> -c ${PGO_CONFIG_FILE} -o ${CMAKE_BINARY_DIR}/pgo_run -r
    DEPENDS ${PROJECT_NAME}
    COMMENT "Running PGO training workload"
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(APPLE)
        execute_process(
            COMMAND xcrun -f llvm-profdata
            OUTPUT_VARIABLE LLVM_PROFDATA_EXECUTABLE
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()
    if(NOT LLVM_PROFDATA_EXECUTABLE)
        find_program(LLVM_PROFDATA_EXECUTABLE llvm-profdata)
    endif()
    if(LLVM_PROFDATA_EXECUTABLE)
        add_custom_target(pgo-merge
            COMMAND ${CMAKE_COMMAND} -E make_directory ${PGO_DIR}
            COMMAND ${LLVM_PROFDATA_EXECUTABLE} merge
                    -output=${PGO_DIR}/default.profdata
                    ${PGO_DIR}/default.profraw
            DEPENDS pgo-run
            COMMENT "Merging Clang PGO profiles"
        )
    else()
        message(STATUS "llvm-profdata not found; Clang PGO merge target disabled.")
    endif()
endif()

# Testing — only build in Debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    file(GLOB TEST_SOURCES "tests/*.cpp" "tests/run/run_tests.cpp")
    add_executable(test_${PROJECT_NAME} ${TEST_SOURCES} $<TARGET_OBJECTS:shared_objs>)
    target_link_libraries(test_${PROJECT_NAME} PRIVATE ${COMMON_LIBS})
    target_include_directories(test_${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/tests/run)
endif()
