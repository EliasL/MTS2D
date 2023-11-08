# Define the build directory
BUILD_DIR = build

# Target executable file
TARGET = $(BUILD_DIR)/main


ALGLIB_DIR = libs/alglib-cpp
CXXFLAGS += -O3 # Optimization value. Use -O3 for high performance
CXXFLAGS += -ggdb -Xpreprocessor -fopenmp -std=c++14 -fPIC 
CXXFLAGS += -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
LIBFLAGS += -I/opt/homebrew/opt/libomp/lib -lomp
LIBFLAGS += -L/opt/homebrew/opt/llvm/lib/c++ -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++

CPPFLAGS += -isystem $(ALGLIB_DIR)/src #-isystem is the same as -I except it ignores warnings

# Get a list of all subdirectories under src/
SUBDIRS := $(shell find src/ -type d)

# Prefix each subdirectory with -I to add to the include path
INCLUDES := $(addprefix -I,$(SUBDIRS))

# Now add INCLUDES to your CPPFLAGS
CPPFLAGS += $(INCLUDES)

CXX = xcrun  clang++

# This includes all cpp files in the libraries and src folder
SRC += $(wildcard $(ALGLIB_DIR)/src/*.cpp)
SRC += $(shell find src -name '*.cpp') # Cpp and tpp



OBJ = $(SRC:.cpp=.o)

# Rule to compile .cpp files to .o files and place them in the build directory
$(BUILD_DIR)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)  # Create the build directory if it doesn't exist
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Update the dependency on $(BUILD_DIR)/%.o in the target rule
$(TARGET): $(OBJ:%=$(BUILD_DIR)/%)
	@echo Linking $@
	@echo $(CPPFLAGS)
	@mkdir -p $(@D)  # Create the build directory if it doesn't exist
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJ:%=$(BUILD_DIR)/%) $(LIBFLAGS) $(LIBGFLAGS) -o $@
	@echo Build Complete

clean:
	@if [ -d "$(BUILD_DIR)/src" ]; then \
		find $(BUILD_DIR)/src -name '*.o' -delete; \
		rm -f $(BUILD_DIR)/main; \
	fi
	@echo My object files and the main binary removed. Library files are preserved.

cleanAll:
	@rm -rf $(BUILD_DIR)
	@echo All object files and binaries removed.