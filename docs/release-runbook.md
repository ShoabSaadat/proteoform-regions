# proteoform-regions — Release Runbook (Shoaib's click list)

> Generated 2026-10-07 after P7-prep. Everything is verified and staged:
> PR #1 open & mergeable · CI green (run 37609577116) · tag `v0.1.0` local-only ·
> wheel+sdist built · `twine check` passed · docs site built · repo PRIVATE.
> **Nothing below has been executed.** Do the steps in order; each has a verification.

Estimated total time: ~20 minutes (plus wait for TestPyPI/PyPI propagation).

---

## Step 0 — Prerequisites (one-time accounts)

- [ ] Logged into GitHub as **ShoabSaadat**
- [ ] Account on **test.pypi.org** (separate account allowed) AND **pypi.org** —
      same username `ShoabSaadat` on both keeps Step 1 simple

## Step 1 — PyPI trusted publishing (2 pending publishers)

PyPI no longer wants your passwords/API tokens for GitHub Actions — you register
the *workflow* as a trusted publisher instead.

1. Go to `https://test.pypi.org/manage/account/publishing/`
2. Click **Add a new pending publisher** and fill EXACTLY:
   - PyPI project name: `proteoform-regions`
   - Owner: `ShoabSaadat`
   - Repository: `proteoform-regions`
   - Workflow name: `publish.yml`
   - Environment name: `testpypi`
3. Save. Repeat 1–3 on `https://pypi.org/manage/account/publishing/` with
   Environment name: `pypi`

**Verify:** both sites list a pending publisher for `proteoform-regions`.

## Step 2 — GitHub Environments

1. Repo → **Settings → Environments → New environment**
2. Create one named `testpypi` (no protections needed)
3. Create one named `pypi`
   - Recommended: add yourself as required reviewer OR leave open (fine for a
     single-maintainer repo)

**Verify:** Settings → Environments lists `testpypi` and `pypi`.

## Step 3 — TestPyPI dry-run

1. Repo → **Actions** → workflow **publish** → **Run workflow**
2. Branch: `feat/docs-walkthrough-release` (the PR branch — it contains
   `publish.yml`) → **Run workflow**
3. Wait for the job (~2 min).

**Verify:** job green; then check
`https://test.pypi.org/project/proteoform-regions/` exists and
`pip install -i https://test.pypi.org/simple/ --no-deps proteoform-regions==0.1.0`
works in a scratch venv.

> If the job fails on "trusted publisher" validation: re-check Step 1 fields
> character-for-character (owner/repo/workflow/environment are case-sensitive).

## Step 4 — Gate-P7 review (Hermes, before the real release)

Say the word and Hermes will: inspect the TestPyPI artifact, re-run the
installed-CLI smoke test from the TestPyPI wheel, confirm README renders on the
TestPyPI page, and green-light the real release.

## Step 5 — The real release

1. Repo → **Pull requests** → **PR #1** ("P6+P7-prep: fragpipe adapter, Sheet12
   walkthrough, MkDocs site, release workflows") → **Merge** (merge commit is fine)
2. From your machine:
   ```bash
   cd ~/git-repos/proteoform-regions
   git checkout main && git pull
   git push origin v0.1.0        # pushes the prepared tag -> triggers publish.yml
   ```
3. Watch **Actions → publish** on `main`.

**Verify:** job green; `https://pypi.org/project/proteoform-regions/` live;
in a scratch venv: `pip install proteoform-regions` then
`proteoform-regions --version` prints `0.1.0`.

## Step 6 — Docs site (GitHub Pages)

1. Repo → **Settings → Pages** → Source: **GitHub Actions**
2. Repo → **Actions** → workflow **deploy-pages** → **Run workflow** (branch `main`)

**Verify:** `https://shoaabsaadat.github.io/proteoform-regions/` loads with the
MkDocs Material site (may take ~1 min).

## Step 7 — Zenodo DOI (versioned + concept DOI)

1. `https://zenodo.org` → account → **GitHub** section → enable sync for
   `ShoabSaadat/proteoform-regions` (flip the toggle ON)
2. Repo → **Releases → Draft a new release** → choose tag `v0.1.0` →
   title `v0.1.0 — first public release` → any short notes → **Publish release**
3. Zenodo mints a version DOI automatically (first release also creates the
   concept DOI). Grab both from zenodo.org → your record.

**Verify:** the Zenodo record shows the repo, authors, and `v0.1.0` files.

## Step 8 — Repo visibility (the gate to NM submission)

Repo → **Settings → General → Danger Zone → Change visibility → Public**.

> Code availability requires a public repo at submission time. Do this at (or
> just before) the Nature Methods submission, not before — your call entirely.

## Step 9 — Back to the paper (Hermes does these with your go)

- [ ] Insert `https://github.com/ShoabSaadat/proteoform-regions` + PyPI name +
      Zenodo DOI into the manuscript **Code availability** section (currently
      holds "[repository URL to be inserted]")
- [ ] Update **SUBMISSION_README.md** + **submission_package** checklists
- [ ] Add the package to the Barnes response letter (point 4: "done, here's the DOI")
- [ ] Update the explainer artifact's software section

---

## Rollback notes

- TestPyPI release can simply be abandoned (no consequences).
- A bad real release: DO NOT delete on PyPI (versions are permanent); yank is
  available per-version in project settings; next version is 0.1.1.
- Pages/Zenodo/visibility are all reversible independently.
