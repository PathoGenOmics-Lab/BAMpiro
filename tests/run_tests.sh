#!/usr/bin/env bash
#
# Run the whole BAMpiro test suite locally, the same way CI does.
#
#   tests/run_tests.sh            # everything
#   tests/run_tests.sh unit       # only the Python unit tests
#   tests/run_tests.sh js         # only the report front-end tests
#   tests/run_tests.sh pipeline   # only the Nextflow stub runs (needs nextflow + Java 17+)
#   tests/run_tests.sh lint       # only ruff
#   tests/run_tests.sh e2e        # a real run in Docker on a simulated cohort (minutes; not part of 'all')
#
# The Nextflow leg is skipped automatically when the CLI is not installed.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

TARGET="${1:-all}"
FAILED=()

step() { printf '\n\033[1m== %s\033[0m\n' "$1"; }
record() { if "$@"; then :; else FAILED+=("$1"); fi; }

# Nextflow writes a storeDir cache and a results tree; start from a clean slate so the stub
# runs are hermetic and BUILD_MAPPABILITY is actually exercised rather than restored.
rm -rf .test_cache results_test

if [[ "$TARGET" == "all" || "$TARGET" == "lint" ]]; then
    step "ruff"
    if command -v ruff >/dev/null 2>&1; then
        record ruff check .
    else
        echo "ruff not installed, skipping (pip install ruff)"
    fi
fi

if [[ "$TARGET" == "all" || "$TARGET" == "unit" ]]; then
    step "pytest: unit tests for bin/"
    record python3 -m pytest tests/unit -q
fi

if [[ "$TARGET" == "all" || "$TARGET" == "js" ]]; then
    step "node: report front-end"
    if command -v node >/dev/null 2>&1; then
        # Pass the files explicitly: `node --test <dir>` is not supported on every Node version.
        record node --test tests/js/*.test.mjs
    else
        echo "node not installed, skipping"
    fi
fi

if [[ "$TARGET" == "all" || "$TARGET" == "pipeline" ]]; then
    step "pytest: Nextflow stub runs"
    if command -v nextflow >/dev/null 2>&1; then
        record python3 -m pytest tests/pipeline -q
    else
        echo "nextflow not installed, skipping (needs Java 17+)"
    fi
fi

# Not part of 'all': it pulls the pinned images and runs the pipeline for real, which takes minutes.
if [[ "$TARGET" == "e2e" ]]; then
    step "pytest: end-to-end run on the simulated cohort (Docker)"
    if command -v nextflow >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        BAMPIRO_E2E=1 record python3 -m pytest tests/e2e -q -rs
    else
        echo "needs nextflow (Java 17+) and a running Docker, skipping"
    fi
fi

echo
if [[ ${#FAILED[@]} -eq 0 ]]; then
    printf '\033[32mAll checks passed.\033[0m\n'
else
    printf '\033[31mFailed: %s\033[0m\n' "${FAILED[*]}"
    exit 1
fi
