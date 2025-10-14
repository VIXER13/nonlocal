BUILD_DIR := build
TOOLCHAIN_FILE := $(BUILD_DIR)/conan_toolchain.cmake
UNITTEST_FILE := $(BUILD_DIR)/tests/unit_tests
BUILD_MAKEFILE := $(BUILD_DIR)/Makefile
PROFILE_PATH := ./.profiles/clang

# Install conan dependencies if toolchain doesn't exist
$(TOOLCHAIN_FILE): $(PROFILE_PATH)
	conan install . --build=missing --output-folder=$(BUILD_DIR) --profile $(PROFILE_PATH)

# Generate CMake build files if Makefile doesn't exist
$(BUILD_MAKEFILE): $(TOOLCHAIN_FILE)
	cd $(BUILD_DIR) && cmake .. --preset conan-release

# Default target - ensure everything is set up
.PHONY: setup
setup: $(BUILD_MAKEFILE)

# Build target
.PHONY: build
build:
	cmake --build $(BUILD_DIR) --config Release --parallel -- -s

.PHONY: run-tests
run-tests:
	cd $(BUILD_DIR) && ctest --output-on-failure

# Clean target
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

# Help target
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  setup       - Create build directory and configure project"
	@echo "  build       - Build the project"
	@echo "  run-tests   - Run unit tests"
	@echo "  clean       - Remove build directory"
	@echo "  help        - Show this help message (default)"

# Default target is build
.DEFAULT_GOAL := help