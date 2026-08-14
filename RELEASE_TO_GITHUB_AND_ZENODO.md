# GitHub and Zenodo release procedure

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

## 2. Create the GitHub repository

Create an empty repository under the confirmed owner in the GitHub web
interface. Do not ask GitHub to add a README, `.gitignore`, or license because
these are already handled locally. The public repository must be public before
Zenodo can archive it.

Then add the returned HTTPS remote and push the committed `main` branch:

```powershell
git remote add origin https://github.com/OWNER/pdac-spatial-program-association.git
git push -u origin main
```

Replace `OWNER` only with the confirmed GitHub account or organization name.

## 3. Create the tagged GitHub release

Update `CITATION.cff`, `zenodo.json`, and `LICENSE` from their templates only
after author and license metadata are final. Commit those edits, then create a
versioned annotated tag and push it:

```powershell
git add CITATION.cff zenodo.json LICENSE release_manifest.json
git commit -m "Prepare public release v1.0.0"
git tag -a v1.0.0 -m "PDAC spatial program association v1.0.0"
git push origin main --tags
```

Create a GitHub Release from the `v1.0.0` tag. Attach no raw GEO data.

## 4. Archive with Zenodo and obtain the DOI

Sign in to [Zenodo](https://zenodo.org/), open the GitHub integration, enable
the public repository, and select the `v1.0.0` release. Zenodo will archive
the tagged snapshot and mint a version-specific DOI. Check that its title,
authors, license, description, keywords, repository URL, and version match the
final metadata before publishing the record.

Only after Zenodo reports the DOI should the authors add it to the manuscript,
README, and citation metadata in a follow-up tagged commit. Preserve both the
version-specific DOI and Zenodo concept DOI if supplied.

## 5. Post-release integrity check

Confirm that GitHub shows the `v1.0.0` tag, Zenodo resolves the DOI, GitHub
Actions passed static validation, and the archived release contains no raw GEO
archives. Save the GitHub release URL and Zenodo record URL in the submission
cover letter or Data Availability statement.
