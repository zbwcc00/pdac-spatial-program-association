# GitHub and Zenodo v1.0.6 integrity-correction release procedure

## 1. Final local check

From the repository root, activate the canonical Conda environment and run:

```powershell
python scripts/build_release_manifest.py
python scripts/verify_release.py
git status
```

For a complete recomputation with data present, use `RUN_ALL.ps1` as described
in the README, then repeat the two verification commands above.

Review `git status` before staging. It must not contain `data/00_raw`,
`data/01_unpacked`, a local Conda environment, credentials, or `.part` files.

## 2. Update the existing public GitHub repository

The public repository already exists at
`https://github.com/zbwcc00/pdac-spatial-program-association`. Commit and push
only the reviewed v1.0.6 revision files. Do not add raw GEO archives, local
temporary files, or `qa_*` rendering artifacts. The repository must remain
public for Zenodo to archive the release.

Before publishing, confirm that the GitHub file list contains the current
manuscript, final supplement, public runbook, `release_manifest.json`, and
`release_verification.json`, and that all version labels read `v1.0.6`.

```powershell
git push -u origin main
```

## 3. Create the tagged GitHub release

Update `CITATION.cff`, `zenodo.json`, `README.md`,
`release_manifest.json`, and `release_verification.json` with the reviewed
v1.0.6 content. Commit those edits, then create a versioned annotated tag and
push it:

```powershell
git add -A
git commit -m "Release v1.0.6 reproducibility and integrity corrections"
git tag -a v1.0.6 -m "PDAC spatial program association v1.0.6"
git push origin main --tags
```

Create a GitHub Release from the `v1.0.6` tag. Attach no raw GEO data.

## 4. Archive with Zenodo and obtain the DOI

Sign in to [Zenodo](https://zenodo.org/), open the GitHub integration, and
confirm that the public repository is enabled. Creating the GitHub `v1.0.6`
release should then create a new Zenodo version automatically. If it does not,
use Zenodo's GitHub integration to archive that release manually. Check that
its title, authors, license, description, keywords, repository URL, and
version match the final metadata before publishing the record.

Only after Zenodo reports the DOI should the authors add it to the manuscript,
README, and citation metadata in a follow-up tagged commit. Preserve both the
version-specific DOI and Zenodo concept DOI if supplied.

## 5. Post-release integrity check

Confirm that GitHub shows the `v1.0.6` tag, Zenodo resolves the DOI, GitHub
Actions passed static validation, and the archived release contains no raw GEO
archives. Save the GitHub release URL and Zenodo record URL in the submission
cover letter or Data Availability statement.
