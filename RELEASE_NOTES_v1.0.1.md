# Version 1.0.1 — release-integrity correction

This corrective release makes no change to the analytical code, derived
results, figures, manuscript content, input registry, or scientific claims.

It fixes the public integrity manifest so that it can be verified after a
GitHub/Zenodo source archive is unpacked on another platform:

- exclude `__pycache__` directories at every nesting level;
- hash only the LF-normalized, source-controlled release files; and
- fail verification when a versioned text file contains CRLF newlines.

The preceding version-specific Zenodo DOI is
[`10.5281/zenodo.21951888`](https://doi.org/10.5281/zenodo.21951888).
