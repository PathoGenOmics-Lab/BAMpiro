.PHONY: docs docs-serve

# Build the documentation site (HTML) into ./site with MkDocs Material.
# The Markdown lives in docs/*.md (readable from GitHub); config is mkdocs.yml.
# One-off setup: pip install -r docs/requirements.txt
docs:
	mkdocs build --strict

# Build + serve the docs locally with live reload (http://127.0.0.1:8000).
docs-serve:
	mkdocs serve
