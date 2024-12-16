# Makefile

# Compiler and flags for macOS
CXX = clang++
CXXFLAGS = -std=c++17 -Wall -I./include -O3 -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
LDFLAGS = -L/opt/homebrew/opt/libomp/lib -lomp

# Directories
SRC_DIR = src
OBJ_DIR = obj/src
TEST_DIR = test
INCLUDE_DIR = include

# Source files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

# Test files
TEST_TARGETS = ebc ebc_tp_multi_simu

# Default target
all: $(TEST_TARGETS)

# Linking the test executables
ebc: $(OBJ_FILES) $(TEST_DIR)/ebc.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

ebc_tp_multi_simu: $(OBJ_FILES) $(TEST_DIR)/ebc_tp_multi_simu.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# Compiling source files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INCLUDE_DIR)/%.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Generate dependencies
$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	@$(CXX) $(CXXFLAGS) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

# Include dependencies
-include $(OBJ_FILES:.o=.d)

# Clean target
clean:
	rm -rf $(OBJ_DIR)/* $(TEST_TARGETS)

# Debug target to print variables
debug:
	@echo "Source files: $(SRC_FILES)"
	@echo "Object files: $(OBJ_FILES)"
	@echo "Test targets: $(TEST_TARGETS)"

.PHONY: all clean debug
