ROOT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

.PHONY: build
build: docs ## build all binaries
	@mkdir -p $(ROOT_DIR)/cmake-build-debug
	@cd $(ROOT_DIR)/cmake-build-debug && cmake -DCMAKE_BUILD_TYPE=Debug $(ROOT_DIR)
	@cmake --build $(ROOT_DIR)/cmake-build-debug
	@mkdir -p $(ROOT_DIR)/cmake-build-release
	@cd $(ROOT_DIR)/cmake-build-release && cmake -DCMAKE_BUILD_TYPE=Release $(ROOT_DIR)
	@cmake --build $(ROOT_DIR)/cmake-build-release

.PHONY: docs
docs: ## generate full documentation
	@doxygen $(ROOT_DIR)/docs/Doxyfile

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
default: help
