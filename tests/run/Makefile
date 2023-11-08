CXX := xcrun clang++
CXXFLAGS := -ggdb -std=c++14 

# Define the build directory from the main Makefile
BUILD_DIR = ../../build

# List of individual test files (excluding the main runner)
TEST_FILES := $(wildcard ../*_test.cpp)
TEST_OBJECTS := $(patsubst ../%.cpp, %.o, $(TEST_FILES))

# Executable name for your tests
TEST_EXECUTABLE := run_tests

# Main runner source file
MAIN_TESTS := run_tests.o

# The object files from your main project
MAIN_OBJ := $(shell find $(BUILD_DIR)/libs -name '*.o')
MAIN_OBJ += $(shell find $(BUILD_DIR)/src -name '*.o')
MAIN_OBJ := $(filter-out %/main.o, $(MAIN_OBJ))

# Build tests without running them
build_without_running: $(TEST_EXECUTABLE)
	@echo "Tests have been built."

# Compile and run tests
all: $(TEST_EXECUTABLE)
	@./$(TEST_EXECUTABLE) --no-intro=true --no-version=true --no-exitcode=true

# Compile individual test object files
%.o: ../%.cpp
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Link all test objects with the main runner and the already compiled object files
$(TEST_EXECUTABLE): $(MAIN_TESTS) $(TEST_OBJECTS) 
	@$(CXX) $(CXXFLAGS) $^ $(MAIN_OBJ) -o $@ 

# Clean generated files
clean:
	@rm -f $(TEST_EXECUTABLE) $(MAIN_TESTS) $(TEST_OBJECTS)

.PHONY: all clean
