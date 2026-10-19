# Repository Guidelines

## Layout

- `irispy/` is the import package; the distribution name is `irispy-lmsal`.
- Main public data APIs live in `irispy/sji.py` and `irispy/spectrograph.py`.
- Level 2 FITS readers live in `irispy/io/`; analysis and calibration helpers live in `irispy/utils/`.
- Packaged calibration and test data live in `irispy/data/`.
- Tests stay beside their code in `irispy/tests/`, `irispy/io/tests/`, and `irispy/utils/tests/`.
- Documentation is in `docs/`, runnable examples are in `examples/`, and Towncrier fragments are in `changelog/`.
- Do not hand-edit generated files such as `irispy/_version.py`, `docs/api/`, or `docs/generated/`.

## Development Commands

The project supports Python 3.13-3.14; Python 3.14 matches the primary CI job.

```bash
micromamba create -n irispy python=3.14
micromamba activate irispy
python -m pip install -e ".[dev]"
```

Run the full offline suite, including documentation doctests:

```bash
pytest
```

Run network and figure tests explicitly:

```bash
pytest --remote-data=any
pytest --pyargs irispy -m "mpl_image_compare" --mpl --remote-data=any
```

The figure hash library is made on Linux by CircleCI; on macOS the figure tests fail against it, which is not a bug.

CI-style checks are `tox -e py314`, `tox -e codestyle`, and `tox -e build_docs`. To remove generated documentation,
run `make clean` from `docs/`.

## Style and Tests

- Target Python 3.13 syntax. Ruff formats with 120-column lines, double quotes, spaces, and NumPy docstrings.
- Let the `isort` pre-commit hook order imports; Ruff's import-sorting rules are disabled.
- Run `pre-commit run --all-files` before opening a PR.
- Warnings are errors. The default suite excludes `mpl_image_compare` tests and enables RST doctests.
- Mark data downloads with `@pytest.mark.remote_data` and other connectivity-dependent tests with
  `@pytest.mark.online`.
- Add the smallest focused regression test near changed code. Use the existing `figure_test` helper for visual output.

## Scope and Judgment

- Handle only real Level 2 data. No guards, fixtures or code paths for files that do not occur in the archive.
- Slicing and indexing behave like plain `ndcube`.
- Test behavior users see. No tests of import mechanics or other internals, and no xfail tests for workflows
  irispy does not support.
- Pickling a cube is not a supported workflow. Pass arrays across process boundaries; ASDF is the planned save
  format.
- Documentation, gallery text and docstrings say what irispy does and what the user needs. How irispy works
  something out, where a rule came from, why a dependency exists and what is planned belong in code comments or
  the pull request.
- No function that only wraps astropy, ndcube or sunpy; write a recipe in the docs instead. Quick-look and GUI
  tools belong in glue-solar.
- Scientific constants, thresholds and line lists need a citation or a data file.
- Releases are made by people. Agents do not tag, publish or edit release metadata.

## Commits and Pull Requests

- Describe user-visible impact in the PR and link the relevant issue or PR when one exists.
- For user-visible changes, add a ReST Towncrier fragment named `changelog/<PR_NUMBER>.<TYPE>.rst`.
- Valid fragment types are `breaking`, `deprecation`, `removal`, `feature`, `bugfix`, `doc`, and `trivial`.
- Add a counter for multiple fragments of the same type, for example `123.feature.1.rst`.
