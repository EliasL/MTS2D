# Define the build directory
BUILD_DIR := build
ALGLIB_DIR := libs/alglib-cpp

# Target executable file
TARGET := $(BUILD_DIR)/main

# Compiler and Linker
CXX := xcrun clang++

# Compiler flags
CPPFLAGS := -isystem $(ALGLIB_DIR)/src
CXXFLAGS := -O3 -ggdb -Xpreprocessor -fopenmp -std=c++14 -fPIC \
            -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk \
            -MMD -MP
LDFLAGS  := -L/opt/homebrew/opt/llvm/lib/c++ -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++
LDFLAGS  += -L/opt/homebrew/opt/libomp/lib 
LDLIBS := -lomp

# Include directories
SUBDIRS := $(shell find src/ -type d)
INCLUDES := $(addprefix -I,$(subst //,/, $(SUBDIRS)))
CXXFLAGS += $(INCLUDES)

# Sources and objects
SRC := $(wildcard $(ALGLIB_DIR)/src/*.cpp) $(shell find src -name '*.cpp')
SRC := $(subst //,/, $(SRC))
OBJ := $(SRC:%.cpp=$(BUILD_DIR)/%.o)
DEP := $(OBJ:.o=.d)

# Include the dependency files
-include $(DEP)

# Rule to create object files
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Rule to link the final executable
$(TARGET): $(OBJ)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@
	@echo "Build Complete"

clean:
	@find $(BUILD_DIR)/src -type f -name '*.o' -delete
	@find $(BUILD_DIR)/src -type f -name '*.d' -delete
	@$(RM) $(TARGET)
	@echo "Source objects removed."

cleanAll:
	@$(RM) -r $(BUILD_DIR)
	@echo "Source and library objects removed."

all: $(OBJ)

.PHONY: all clean cleanAll
