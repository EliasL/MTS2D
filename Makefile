# Define the build directory
BUILD_DIR = build

# Target executable file
TARGET = $(BUILD_DIR)/gael.exe


ALGLIB_DIR = libs/alglib-cpp

CXXFLAGS += -O3 -ggdb -Xpreprocessor -fopenmp -std=c++14 -fPIC -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk
LIBFLAGS += -I/opt/homebrew/opt/libomp/lib -lomp
LIBFLAGS += -L/opt/homebrew/opt/llvm/lib/c++ -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++

CPPFLAGS += -I$(ALGLIB_DIR)/src

# This includes all cpp files in the libraries and src folder
SRC += $(wildcard $(ALGLIB_DIR)/src/*.cpp)
SRC += $(wildcard src/*.cpp)

OBJ = $(SRC:.cpp=.o)

# Rule to compile .cpp files to .o files and place them in the build directory
$(BUILD_DIR)/%.o: %.cpp
	@echo Compiling $<
	@mkdir -p $(@D)  # Create the build directory if it doesn't exist
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Update the dependency on $(BUILD_DIR)/%.o in the target rule
$(TARGET): $(OBJ:%=$(BUILD_DIR)/%)
	@echo Linking $@
	@mkdir -p $(@D)  # Create the build directory if it doesn't exist
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJ:%=$(BUILD_DIR)/%) $(LIBFLAGS) $(LIBGFLAGS) -o $@
	@echo Build Complete

clean:
	rm -f $(OBJ) $(TARGET)
	@echo All object files and binaries removed
