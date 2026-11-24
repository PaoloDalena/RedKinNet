## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* Possibly misspelled words in DESCRIPTION: "Manichaikul", "et", "al". These are the
  surname of the first author and the abbreviation "et al." in the reference
  Manichaikul et al. (2010) <doi:10.1093/bioinformatics/btq559>; they are spelled correctly.

## Test environments

* local: macOS 26.1 (aarch64), R 4.5.2
* GitHub Actions: ubuntu-latest (R devel, release and oldrel-1), windows-latest (R release),
  macos-latest (R release)
* win-builder: R devel (2026-11-17 r90926) and R release (4.6.2)

## Notes for the reviewers

* The vignettes use only the synthetic example cohort shipped with the package and build in
  about a minute on win-builder.
* Tests that take longer than a few seconds are skipped on CRAN.
