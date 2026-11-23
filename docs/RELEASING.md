# Releasing pypedigree

Publishing is automated by `.github/workflows/release.yml`, which runs on a
version tag. Authentication uses PyPI **Trusted Publishing** (OIDC), so there
is no API token stored anywhere in the repository or in GitHub secrets.

## One-time setup on PyPI

This must be done once, by hand, before the first release. It cannot be
scripted from the repository.

1. Sign in to <https://pypi.org> and go to **Your projects → Publishing**, or
   for a project that does not exist yet:
   <https://pypi.org/manage/account/publishing/>.
2. Add a **pending publisher** (GitHub) with exactly these values:

   | Field | Value |
   |-------|-------|
   | PyPI project name | `pypedigree` |
   | Owner | `bourgeron-lab` |
   | Repository name | `pypedigree` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

   The environment name must match the `environment: name: pypi` in the
   workflow, or PyPI will reject the token. PyPI compares the owner,
   repository and workflow names character for character, case included:
   any difference fails the upload with `invalid-publisher` ("Publisher with
   matching claims was not found").

3. In the GitHub repository, go to **Settings → Environments → New
   environment** and create one named `pypi`. Adding a required reviewer
   there is recommended: it makes every publish pause for an explicit
   approval.

The project name `pypedigree` was free on PyPI at the time of writing. Verify
it still is before setting up the publisher:

```bash
curl -s -o /dev/null -w '%{http_code}\n' https://pypi.org/pypi/pypedigree/json
```

`404` means the name is available; `200` means it is taken.

## Cutting a release

1. Make sure `main` is up to date and CI is green.

2. Move the `[Unreleased]` section of `CHANGELOG.md` under a new version
   heading with today's date:

   ```markdown
   ## [0.5.3] - 2027-01-15
   ```

   This section becomes the GitHub release notes: the release workflow
   extracts it (from the heading to the next `## [`) and fails before
   publishing anything if the version has no section.

3. Bump the version in **one** place — `pyproject.toml`:

   ```toml
   version = "0.5.3"
   ```

   `pypedigree.__version__` reads the installed package metadata, so no
   other code needs editing. Then update the places that print the version:
   `version` and `date-released` in `CITATION.cff`, the two `v0.x.y` strings
   in `docs/pypedigree_documentation.html` (header and footer; regenerate the
   PDF), the placeholder in `.github/ISSUE_TEMPLATE/bug_report.yml`, and the
   project's own entry in `uv.lock` with `uv lock`.

   A release candidate is versioned `1.0.0rc1` (PEP 440) and tagged
   `v1.0.0rc1`, with its own CHANGELOG section `## [1.0.0rc1] - date`. PyPI
   treats it as a pre-release (`pip install pypedigree` keeps the last final
   version; `pip install pypedigree==1.0.0rc1` gets it), the GitHub release is
   marked as a pre-release, and the web app is not republished for it.

   A tag that was pushed is never reused, even if nothing was published:
   bump to the next version instead (0.4.0 and 0.4.1 were withdrawn this
   way and released as 0.4.2).

   Redraw the documentation's images so they show this version:
   `uv run python docs/imgs/make_notation_images.py` and
   `uv run python docs/imgs/examples/make_example_images.py`.

4. Verify locally:

   ```bash
   uv sync --extra dev --extra pandas
   uv run pytest -q
   uv run ruff check src/ tests/ web/ && uv run ruff format --check src/ tests/ web/
   uv build && uvx twine check --strict dist/*
   ```

5. Commit, tag and push. The tag must be the version prefixed with `v` — the
   workflow fails the build if it disagrees with `pyproject.toml`:

   ```bash
   git commit -am "Release 0.5.3"
   git tag -a v0.5.3 -m "Release 0.5.3"
   git push origin main --follow-tags
   ```

6. The workflow then runs the test suite, builds the sdist and wheel, checks
   the metadata with `twine check --strict`, re-runs the tests from inside the
   extracted sdist, publishes to PyPI, and creates the GitHub release with
   the artifacts attached.

## The web app

The same tag runs `.github/workflows/pages.yml`, which publishes the web app
(`web/`) on GitHub Pages with the wheel of that release:
<https://bourgeron-lab.github.io/pypedigree/>. It needs, once, **Settings →
Pages → Build and deployment → Source: GitHub Actions**. A change to `web/`
alone can be published from **Actions → Web app → Run workflow** on `main`.
The wheel is then built from `main`, so do it only when `src/` has no
unreleased change; otherwise wait for the next release.

## Rehearsing without publishing

Run the workflow manually from the **Actions** tab with `dry_run` left at
`true`. Everything is built and checked, and the publish step is skipped.

## Making the release citable

The repository is connected to
[Zenodo](https://zenodo.org/account/settings/github/): each published GitHub
release is archived there and gets its own version DOI, with no action needed.
`CITATION.cff` and the README cite the concept DOI,
[10.5281/zenodo.24289471](https://doi.org/10.5281/zenodo.24289471), which covers every version and
resolves to the latest, so they need no change per release. 0.4.4 was the
first release archived (10.5281/zenodo.24289472).

If a release does not appear on Zenodo, check the webhook deliveries in
**Settings → Webhooks**: Zenodo acts on the `release` event with action
`published` (a 500 on the `released` action is harmless).

GitHub reads `CITATION.cff` and shows a **Cite this repository** button on the
repository page.
