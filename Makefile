COMPILER ?= clang
BUILD_DIR := build
TOOLCHAIN_FILE := $(BUILD_DIR)/conan_toolchain.cmake
UNITTEST_FILE := $(BUILD_DIR)/tests/unit_tests
BUILD_MAKEFILE := $(BUILD_DIR)/Makefile
PROFILE_PATH := ./.profiles/$(COMPILER)
COMPILER_MARKER := $(BUILD_DIR)/.compiler.stamp

$(COMPILER_MARKER):
	mkdir -p $(BUILD_DIR)
	echo $(COMPILER) > $@

$(TOOLCHAIN_FILE): $(COMPILER_MARKER) $(PROFILE_PATH)
	conan install . --build=missing --output-folder=$(BUILD_DIR) --profile $(PROFILE_PATH)

$(BUILD_MAKEFILE): $(COMPILER_MARKER) $(TOOLCHAIN_FILE)
	cd $(BUILD_DIR) && cmake .. --preset conan-release

.PHONY: update_compiler
update_compiler:
	@[ "$$(cat $(COMPILER_MARKER))" != "$(COMPILER)" ] && echo "$(COMPILER)" > $(COMPILER_MARKER) || true

# Setup target
.PHONY: setup
setup: $(BUILD_MAKEFILE)

# Build target
.PHONY: build
build: update_compiler setup
	cmake --build $(BUILD_DIR) --config Release --parallel -- -s

# Run unit tests
.PHONY: run-tests
run-tests:
	cd $(BUILD_DIR) && ctest --output-on-failure

# Clean target
.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

# Print help
.PHONY: help
help:
	@echo "Available targets:" 
	@echo "  setup         - Configure project (reconfigures if COMPILER changed)" 
	@echo "  build         - Build the project" 
	@echo "  run-tests     - Run unit tests (builds first)" 
	@echo "  clean         - Remove build directory" 
	@echo "  help          - Show this help message (default)" 

.DEFAULT_GOAL := help