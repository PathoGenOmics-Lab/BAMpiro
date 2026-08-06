## What changed and why

<!-- A couple of sentences. Link the issue if there is one. If this changes a default or an output
     file, say so explicitly: that is what reviewers most need to notice. -->

## How it was verified

<!-- Paste the result, or say which legs you ran. Note anything you could not run locally, e.g. the
     stub run without Nextflow installed. -->

```bash
tests/run_tests.sh
```

## Checklist

- [ ] `tests/run_tests.sh` passes locally (or the legs I could run do, and I said which)
- [ ] `ruff check .` is clean, and I did not restyle code the change does not touch
- [ ] New or changed behaviour has a test, and tests pin what the code does today
- [ ] A new process in `modules/` has a `stub:` block, and `-profile test_full` reaches it
- [ ] `tests/data/` was not edited by hand (fixtures come from `make_test_data.py`)
- [ ] User-visible changes have a `CHANGELOG.md` entry and, where relevant, a docs update
