.PHONY: lint install-hooks

# Run the project's pinned Ruff configuration from the local Python environment.
lint:
	cd python && .venv/bin/ruff check tskit_extras/ tests/

# GitHub Desktop uses this repository's Git hooks when committing.
install-hooks:
	git config core.hooksPath .githooks
