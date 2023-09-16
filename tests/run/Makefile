CXX := clang++
CXXFLAGS := -std=c++14

# List of individual test files (excluding the main runner)
TEST_FILES := $(wildcard ../*_test.cpp)
TEST_OBJECTS := $(patsubst tests/%.cpp, tests/%.o, $(TEST_FILES))

# Executable name for your tests
TEST_EXECUTABLE := run_tests

# Main runner source file
MAIN_TESTS := run_tests.cpp

# Compile and run tests
all: $(TEST_EXECUTABLE)
	./$(TEST_EXECUTABLE)

# Compile individual test object files
tests/%.o: tests/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link all test objects with the main runner
$(TEST_EXECUTABLE): $(MAIN_TESTS) $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ 

# Clean generated files
clean:
	rm -f $(TEST_EXECUTABLE) $(TEST_OBJECTS)

.PHONY: all clean
