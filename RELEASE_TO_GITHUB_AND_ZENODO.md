# GitHub and Zenodo v1.0.8 document-consistency release procedure

## Purpose

V1.0.8 corrects only version-specific archive references in the main manuscript
and supplementary materials. Its reserved Zenodo version DOI is
`10.5281/zenodo.22032109`. It does not change analytical outputs, figures,
input-data scope, or scientific conclusions.

## Local integrity gate

From the repository root, run:

```powershell
python scripts/build_release_manifest.py
python scripts/verify_release.py
```

The verification report must show `status: PASS` with zero failures. Do not
stage raw GEO archives, local environments, caches, `.part` files, or QA
artifacts.

## Git tag

Commit the reviewed files, then create and push an annotated tag:

```powershell
git add -A
git commit -m "Release v1.0.8 document-version consistency correction"
git tag -a v1.0.8 -m "PDAC spatial program association v1.0.8"
git push origin main
git push origin v1.0.8
```

## Zenodo draft

The v1.0.8 Zenodo draft was created manually to reserve its version DOI before
document assembly. Build the source-and-results ZIP from a freshly extracted
`v1.0.8` Git archive, verify it, and upload that one verified ZIP to the
existing draft. Do not import v1.0.7 files, and do not upload raw GEO data.

Before publishing, confirm that the record version is `1.0.8`, the DOI is
`10.5281/zenodo.22032109`, visibility is public, and the record metadata match
`zenodo.json`. Publish the existing draft only after checking the uploaded ZIP.

Because this Zenodo record is a manual draft, do not create a separate GitHub
Release until confirming that the GitHub-Zenodo integration will not create a
duplicate record.
