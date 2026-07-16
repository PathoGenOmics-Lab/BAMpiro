.PHONY: docs docs-serve

# Build the documentation site (HTML) into docs/book/ with mdbook.
# The GitHub-facing Markdown lives in docs/*.md; the mdbook sources are in docs/src/.
docs:
	cd docs && mdbook build

# Build + serve the docs locally with live reload (http://localhost:3000).
docs-serve:
	cd docs && mdbook serve --open
