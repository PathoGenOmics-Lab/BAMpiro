# Contributing to BAMpiro

BAMpiro is a Nextflow (DSL2) pipeline maintained by a small group at I2SysBio (University of
Valencia - CSIC) as part of the [PathoGenOmics-Lab](https://github.com/PathoGenOmics-Lab). Bug
reports, questions and patches are all welcome, and none of them need to be large: a corrected
default, a clearer error message or a missing line in the docs is a useful contribution.

This document covers how to get a working setup, how to run the tests, what CI will check, and the
few local conventions that are not obvious from reading the code.

Everything here happens under the [Code of Conduct](CODE_OF_CONDUCT.md). Contributions are accepted
under the project's [GPL-3.0 license](LICENSE).

## Getting set up

You need two things to run the pipeline itself:

- **Nextflow** `>= 24.04.2`, which needs **Java 17+**
- **Docker** or **Singularity / Apptainer** for the container

The cleanest way to get both without touching your system Java is a conda environment, which is what
[Your first run](docs/tutorials/first-run.md) recommends:

```bash
conda create -n bampiro -c conda-forge -c bioconda nextflow
conda activate bampiro
nextflow -version   # confirm it is >= 24.04.2
```

Then clone the repository and check that it runs:

```bash
git clone https://github.com/PathoGenOmics-Lab/BAMpiro
cd BAMpiro
nextflow run main.nf -profile test -stub-run
```

That last command executes the whole DAG against the bundled 170 kB fixture cohort in a couple of
seconds. It pulls no container, reads no reference genome and touches the network only if Nextflow
needs to fetch a plugin.

For the test suite and the linter you also want:

- **Python** `>= 3.10` with `pytest` and `numpy` (`pip install pytest numpy ruff`)
- **Node** `>= 18` for the report front-end tests

Working on the docs additionally needs MkDocs Material: `pip install -r docs/requirements.txt`, then
`make docs` (a strict build) or `make docs-serve` (live reload on http://127.0.0.1:8000).

## Running the tests

One script runs everything, the same way CI does:

```bash
tests/run_tests.sh              # everything
tests/run_tests.sh lint         # ruff only
tests/run_tests.sh unit         # Python unit tests for bin/ only
tests/run_tests.sh js           # report front-end only
tests/run_tests.sh pipeline     # Nextflow stub runs only
```

There are around 1,100 tests in three legs:

| Leg | What it covers | Needs |
| :--- | :--- | :--- |
| `tests/unit/` | The Python under `bin/`: the consensus decision tree, the QC verdict engine, the parsers, the k-mer liftover, the read filter | `pytest`, `numpy` |
| `tests/js/` | The report's hand-written ES5 statistics, checked against SciPy and statsmodels reference values, plus the integrity of the asset bundle | `node >= 18` |
| `tests/pipeline/` | Samplesheet validation and a full `-stub-run` of the DAG on both test profiles | `nextflow`, Java 17+ |

A leg whose tool is not installed skips itself instead of failing, so a partial local setup still
gives useful signal. Nothing needs a container, a reference genome, real reads or a network
connection. [`tests/README.md`](tests/README.md) documents the fixture cohort, the two test profiles
and the Node harness in more detail.

Two things to know before you add a test:

- Anything written to disk goes in `tmp_path`. Do not add files to `tests/data/` by hand. If a
  fixture is worth keeping, add it to `tests/data/make_test_data.py`, which derives the whole cohort
  from one seed. CI regenerates the cohort and fails on any diff, so a hand-added file breaks the
  build.
- Unit tests import a script from `bin/` with `load_script("name")` from `tests/conftest.py`; `bin/`
  is not a package. Importing `qc_report` reads `bin/report_assets/`, so it only works in place.

## What CI runs

[`.github/workflows/ci.yml`](.github/workflows/ci.yml) gates every pull request with five jobs:

| Job | What it does |
| :--- | :--- |
| Lint (ruff) | `ruff check .` |
| Unit tests | `pytest tests/unit` on Python 3.10 and 3.13 |
| Fixtures are reproducible | Regenerates `tests/data/` and fails on any diff |
| Report front-end (Node) | `node --test tests/js/*.test.mjs` on Node 20 |
| Stub run | `pytest tests/pipeline` against Nextflow `24.04.2` and `latest-stable` |

The two Python versions are the oldest the scripts must still run on and the one in the container.
The two Nextflow versions are the minimum the manifest declares and the current release, the latter
as an early warning rather than a hard requirement.

[`.github/workflows/docs.yml`](.github/workflows/docs.yml) additionally runs `mkdocs build --strict`
on pull requests that touch `docs/`, `mkdocs.yml` or the workflow itself, so a broken link or a page
missing from the nav fails review rather than `main`.

Actions minutes are metered on this repository, so a new push cancels the previous run of the same
branch. Run `tests/run_tests.sh` locally before you push rather than using CI as your test runner.

## Linting: narrow on purpose

`ruff` is configured in [`pyproject.toml`](pyproject.toml) with a deliberately narrow rule set:

```toml
select = ["F", "E9", "B006", "B008", "B017", "B904"]
line-length = 120
```

That catches real defects (undefined names, unused imports and variables, broken f-strings, syntax
errors, mutable default arguments) and nothing else. The Python under `bin/` chose a compact, dense
layout, and the rule set is narrow so that choice is not relitigated file by file.

**Do not restyle the codebase.** No repository-wide formatter run, no import reordering, no
reflowing of code you are not otherwise changing. A diff that mixes a one-line fix with 400 lines of
reformatting cannot be reviewed. If you think a rule should be added, propose it in an issue on its
own, separately from the change that prompted it.

`bin/extract_kraken_reads.py` is excluded from linting entirely. See
[Two things worth knowing](#two-things-worth-knowing) below.

## Commit messages

The convention in this repository is a lowercase area prefix, a colon, and a short lowercase summary
of what the change does, followed by a blank line and two to four sentences of prose explaining it.
No bullet lists, no headed sections, no trailers. Read `git log --oneline -30` for the shape of it.

The prefix is usually a conventional-commit type (`fix:`, `test:`, `docs:`, `ci:`, `refactor:`), but
the part of the pipeline you touched works just as well and is used often (`report:`, `liftover:`,
`mask:`, `pathotypr:`, `config:`, `container:`). A parenthesised scope is fine where it helps:
`docs(tutorials):`, `ci(docs):`.

The body says why the change is the way it is, not what the diff already shows:

```
ci: validate the docs on pull requests instead of duplicating the build

docs.yml already runs mkdocs build --strict with the system libraries the social
plugin needs, so it just gains a pull_request trigger; only the build job runs
there, since deploy stays gated on the Pages variable. Drops the stale
indel-mask branch trigger now that the branch is merged.
```

## Adding a process

Every process in `modules/` needs a `stub:` block. A stub `touch`es each declared output and nothing
else:

```groovy
stub:
"""
touch ${sampleId}.${refId}.consensus.fasta
touch ${sampleId}.${refId}.consensus.log
"""
```

Without one, `-stub-run` aborts and the whole pipeline leg of CI goes red.
`tests/pipeline/test_stub_run.py::test_there_is_a_stub_block_for_every_process` also checks this
statically, so you get a clear failure instead of a confusing Nextflow error.

If the process sits on an optional branch, make sure `-profile test_full` reaches it. That profile
turns on every optional feature (Kraken, pathotypr typing, canonical annotation, the liftover, the
virgin consensus) precisely so one stub run instantiates every declared process, and
`test_the_full_profile_reaches_every_process` fails when something drops out of that coverage.

## Tests pin observed behaviour

The suite records what the code does **today**, not what it ought to do. When a test looks like it is
asserting something wrong, that is usually deliberate and there is a comment saying so.

If you find a genuine defect, do not leave a failing test behind to mark it. Split it in two:

1. A test that pins the current, wrong behaviour, with a comment explaining why it is wrong.
2. A separate change that fixes the code and flips the test.

Both can be in the same pull request, as separate commits. What matters is that `main` is never red:
a suite that is expected to fail is a suite everyone learns to ignore.

## Changelog

[`CHANGELOG.md`](CHANGELOG.md) follows [Keep a Changelog](https://keepachangelog.com/) and the
project follows [Semantic Versioning](https://semver.org/). Anything a user would notice (a new
parameter, a changed default, a fixed wrong result, a new report panel) gets an entry under
`Added`, `Changed` or `Fixed` in the section for the version being prepared. Add an `## [Unreleased]`
section at the top if there is not one yet.

Internal work with no user-visible effect (refactors, tests, CI, most docs changes) does not need an
entry.

## Two things worth knowing

**`main.nf` has to satisfy the strict Nextflow parser.** From 25.10 that parser is the default, and
it is stricter than the language you may be used to: no statements at the top level of a script
(everything lives in a `workflow` or a function), no `while` loops, no C-style `for`, no assignment
used as an expression, and a dynamic process directive must be a closure -
`publishDir path: { "..." }`, not `publishDir "..."`. If a run dies with
`Statements cannot be mixed with script declarations`, that is what you have hit. CI runs both the
minimum supported 24.04.2 and `latest-stable`, so either parser rejecting your change fails the
build.

**`bin/extract_kraken_reads.py` is vendored and must stay byte-identical to upstream.** It comes
from [KrakenTools](https://github.com/jenniferlu717/KrakenTools) (the URL is `params.krakentools_url`
in `nextflow.config`, and `modules/qc.nf` falls back to downloading it when the local copy is
missing). It is not ours to restyle, refactor or "fix", and `pyproject.toml` excludes it from ruff
for that reason. Updating it means replacing it with a newer upstream copy, in its own commit.

## Opening a pull request

Branch off `main`, keep the change focused, and fill in
[the pull request template](.github/PULL_REQUEST_TEMPLATE.md). Run `tests/run_tests.sh` before you
push. If you are unsure whether an idea fits, open an issue first and ask; that is cheaper than
writing code that gets turned down.

## Reporting bugs and asking questions

Use the [issue templates](https://github.com/PathoGenOmics-Lab/BAMpiro/issues/new/choose). The bug
template asks for the things that actually make a BAMpiro failure reproducible: your Nextflow
version, the exact command, the profile, the container engine, the failing process name and the
contents of its work directory (`.command.err`, `.command.sh`, `.command.log`). Nextflow prints the
work directory of a failed task, so `cd` there first.

Before you file, check [Troubleshooting](docs/troubleshooting.md). Most first-run failures are one of
the entries in that table.
