# Releasing

A release is a git tag. Everything else follows from it, so the order matters.

## 1. Prepare

- [ ] `tests/run_tests.sh` is green, and CI is green on `main`.
- [ ] `CHANGELOG.md` has a section for the new version with today's date, and the
      `## [Unreleased]` section is empty.
- [ ] The version is bumped in **three** places, which must agree:
      `nextflow.config` (`manifest.version`), `CITATION.cff` (`version`), and the
      version badge in `README.md`.

## 2. Publish the container first

The pipeline pins its image **by digest**, so the image has to exist before the
tag that references it.

```bash
gh workflow run docker-publish.yml -f tag=<version>
```

It builds once and pushes to both GHCR and Docker Hub. Copy the `sha256:` digest
the action prints into `nextflow.config`:

```groovy
container = "docker://paururo/bampiro@sha256:<digest>"
```

Commit that, then verify a run resolves it:

```bash
nextflow run main.nf -profile test -stub-run
```

!!! note
    Whichever registry the cluster pulls from has to be **public**, or Singularity
    needs `APPTAINER_DOCKER_USERNAME` / `APPTAINER_DOCKER_PASSWORD` in the
    environment. Docker Hub is the public one today.

## 3. Tag and release

```bash
git tag <version>
git push --tags
gh release create <version> --title "<version>" --notes-file <notes.md>
```

The tag triggers `docker-publish.yml` again, which is harmless: it rebuilds the
same content and republishes the same tags.

## 4. Archive for citation

Currently **blocked**: Zenodo can only archive a public repository, and this one
is private. When it goes public:

1. Sign in to [Zenodo](https://zenodo.org) with GitHub and enable the switch for
   `PathoGenOmics-Lab/BAMpiro`.
2. Cut a release. Zenodo archives it and mints a DOI, using the metadata in
   `.zenodo.json`.
3. Add the **concept DOI** (the one that always resolves to the newest version,
   not the per-version DOI) to `CITATION.cff` as `doi:` and to the README badge
   row.

## 5. Docs

`docs.yml` publishes GitHub Pages on push to `main`, gated on the repository
variable `ENABLE_PAGES`. Pages on a private repository needs a paid plan, so
while this repository is private the site is off and the README's docs badges do
not resolve. Nothing to do at release time beyond checking the build passed.
