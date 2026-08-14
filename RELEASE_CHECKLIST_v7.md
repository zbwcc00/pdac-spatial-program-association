# V7 local release checklist

- [x] The spatial primary endpoint, conditional estimand, and randomization hierarchy are recorded.
- [x] Global-Moran and mask-Moran graph procedures are labelled as complementary sensitivities.
- [x] GSE277116 is labelled as sample-level technical replication only.
- [x] GSE202051 is labelled as author-annotated program-attribution/specificity QC only.
- [x] GSE202051 two-sided patient-level test outputs, fixed program lists, control genes, and archive hash are retained.
- [x] GSE202051 has been added to Supplementary Table S21.
- [x] A local manifest hashes release scripts, key results, manuscript sources, figures, and supplementary assets.
- [x] A verification script checks every listed local hash plus the GSE202051 compressed input hash.
- [ ] Create a public GitHub repository without raw public archives or local temporary files.
- [ ] Add an environment lock or package export appropriate to the final compute environment.
- [ ] Tag the exact reviewed commit.
- [ ] Archive that tag with Zenodo and obtain a DOI.
- [ ] Replace local-release placeholders in the manuscript only after the public archive is verified.
