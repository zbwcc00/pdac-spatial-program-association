# v1.0.6 integrity-correction release

## Purpose

This release corrects reproducibility packaging defects identified in the
v1.0.5 archive. It does not alter the study's biological interpretation or
primary estimand.

## Corrections

- The public release manifest excludes Python bytecode and `__pycache__`
  directories at any depth.
- The release verifier relies on exact recorded byte counts and SHA-256 hashes;
  it does not impose a separate line-ending rule.
- Seeded array-block surrogates traverse stable, sorted block identifiers,
  eliminating hash-order-dependent random-number consumption.
- The full run includes the GSE278687 4x4, 8x8, and 12x12 block-size
  sensitivity analysis (999 draws per size; seed 20260813).
- Program scoreability is regenerated from the locked GSE278687/GSE277116
  matrix readers and summed-count scoring function before threshold sensitivity
  analysis; the release no longer relies on an untracked scoreability table.
- Representative GSE277116 spatial maps are regenerated from the locked
  eligible packages using the prespecified within-stratum median-effect rule.

## Release gate

This staging copy is publishable only if `RUN_ALL.ps1` completes and the
newly generated `release_verification.json` reports `"status": "PASS"` with
all required public raw inputs locally verified. `build_public_release_zip.py`
then archives only manifest-listed public files, never raw GEO directories.
