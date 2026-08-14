# Required information before public release

The analysis package is ready to be placed under Git, but the following facts
must be supplied by the corresponding author before GitHub and Zenodo can be
used. They cannot be inferred safely from the local files.

1. Final repository owner and name, for example `github-user/pdac-spatial-program-association`.
2. Public or private GitHub visibility. Zenodo can mint a public DOI only for a public release.
3. Final author names, order, affiliations, corresponding author email, and ORCID iDs where available.
4. A license decision. For openly reusable analysis code, `MIT` or `Apache-2.0` is commonly suitable; this is a legal/author decision, not an automated choice.
5. Final manuscript title, version, and release date.
6. A GitHub login in the browser or an authenticated Git remote. GitHub CLI is not required.
7. A Zenodo login associated with the intended depositor. Either use the Zenodo GitHub integration after pushing the release tag, or provide an authorised personal-access token through a secure local mechanism. Never store a token in this repository.

After these items are confirmed, follow `RELEASE_TO_GITHUB_AND_ZENODO.md`. A DOI must be copied from the completed Zenodo record; it must never be invented or reserved in manuscript text.
