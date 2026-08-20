# v1.0.7 archive-byte-consistency release

## Purpose

This release corrects the byte-level mismatch between the v1.0.6 release
manifest generated in a Windows working tree and the line-ending-normalized
files emitted by a Git source archive. It does not alter any analysis, figure,
manuscript result, primary estimand, or biological interpretation.

## Correction

- The public-file manifest is rebuilt from the exact bytes of the Git archive
  that is released to readers.
- The verifier is run against that archive-derived source tree and against a
  newly built manifest-only ZIP extracted into a clean directory.
- The public ZIP continues to exclude raw GEO inputs and transient Python
  caches.

## Citation boundary

Use the v1.0.7 Git tag and its version-specific Zenodo DOI for final manuscript
submission. The Zenodo concept DOI remains the release-family identifier.
