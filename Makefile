.PHONY: release
release: ## build release binary
	@mkdir -p cmake-build-release
	@cd cmake-build-release && cmake -DCMAKE_BUILD_TYPE=Release ..
	@cd cmake-build-release && cmake --build .
	@./cmake-build-release/irrpolygf2

.PHONY: debug
debug: ## build release binary
	@mkdir -p cmake-build-debug
	@cd cmake-build-debug && cmake -DCMAKE_BUILD_TYPE=Debug ..
	@cd cmake-build-debug && cmake --build .
	@./cmake-build-debug/irrpolygf2

.PHONY: docs
docs: ## generate full documentation
	@cd docs && doxygen Doxyfile

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
default: help
