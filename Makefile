BUILD_DIR := build
TOOLCHAIN_FILE := $(BUILD_DIR)/conan_toolchain.cmake
UNITTEST_FILE := $(BUILD_DIR)/tests/unit_tests
BUILD_MAKEFILE := $(BUILD_DIR)/Makefile
COMPILER_MARKER := $(BUILD_DIR)/.compiler.stamp
THREADS_MARKER := $(BUILD_DIR)/.threads.stamp

DEFAULT_COMPILER := gcc
DEFAULT_THREADS := $(shell command -v nproc > /dev/null 2>&1 && nproc || sysctl -n hw.ncpu)
COMPILER ?= $(if $(wildcard $(COMPILER_MARKER)),$(shell cat $(COMPILER_MARKER)),$(DEFAULT_COMPILER))
THREADS ?= $(if $(wildcard $(THREADS_MARKER)),$(shell cat $(THREADS_MARKER)),$(DEFAULT_THREADS))
PROFILE_PATH := ./.profiles/$(COMPILER)

$(COMPILER_MARKER):
	mkdir -p $(BUILD_DIR)
	echo $(COMPILER) > $@

$(TOOLCHAIN_FILE): $(COMPILER_MARKER) $(PROFILE_PATH)
	conan install . --build=missing --output-folder=$(BUILD_DIR) --profile $(PROFILE_PATH)

$(BUILD_MAKEFILE): $(COMPILER_MARKER) $(TOOLCHAIN_FILE)
	cd $(BUILD_DIR) && cmake .. --preset conan-release

.PHONY: update_compiler
update_compiler:
	@mkdir -p "$(BUILD_DIR)"
	@[ "$$(cat "$(COMPILER_MARKER)" 2>/dev/null)" = "$(COMPILER)" ] || printf '%s\n' "$(COMPILER)" > "$(COMPILER_MARKER)"

.PHONY: update_threads
update_threads:
	@mkdir -p "$(BUILD_DIR)"
	@[ "$$(cat "$(THREADS_MARKER)" 2>/dev/null)" = "$(THREADS)" ] || printf '%s\n' "$(THREADS)" > "$(THREADS_MARKER)"

# Setup target
.PHONY: setup
setup: $(BUILD_MAKEFILE)

# Build target
.PHONY: build
build: update_compiler update_threads setup
	cmake --build $(BUILD_DIR) --config Release -j$(THREADS) -- -s

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