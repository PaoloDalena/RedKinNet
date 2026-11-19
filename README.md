
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RedKinNet <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/PaoloDalena/RedKinNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaoloDalena/RedKinNet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

RedKinNet selects sets of unrelated individuals from pairwise kinship
estimates, such as those of KING or PLINK 2. Pairs related above a
threshold form a network whose connected components, the families, are
resolved exactly, so that the selection is always as large as possible,
also in the large and intricate families of genetic isolates. Among the
many selections of the largest size, RedKinNet can draw one uniformly at
random, with the exact probability that each individual is selected;
choose by preferences for individuals; or reach a target proportion or
mean of a covariate over the whole selection.

## Installation

You can install the development version of RedKinNet from GitHub with:

``` r
# install.packages("remotes")
remotes::install_github("PaoloDalena/RedKinNet")
```

## Example

The package includes `rkn_example`, a synthetic cohort with a kinship
table and covariates.

``` r
library(RedKinNet)

# Check the samples: duplicates and possible contamination
chk <- rkn_check_samples(rkn_example$kinship, threshold = "third",
                         samples = rkn_example$samples$id)

# The network of relatives up to the third degree, without the contaminated sample
g <- rkn_graph_from_kinship(rkn_example$kinship, threshold = "third",
                            samples = rkn_example$samples$id,
                            exclude = chk$contamination$id,
                            covariates = rkn_example$samples, check = FALSE)

# A largest set of unrelated individuals, drawn uniformly among all of them
sel <- rkn_select_random(g, uniform = TRUE, seed = 1)
sel
#> <rkn_selection> Uniform random selection, seed 1
#> Call: rkn_select_random(g = g, uniform = TRUE, seed = 1)
#> Individuals: 751 | selected: 452 (60.2%)
#> Families: 376 | size min/mean/max: 1 / 2.00 / 31
#> Families with up to 5 members: 357 | with 6 or more: 19
#> Families sampled uniformly among their maximum sets: 376 of 376
#> Maximum sets: 4.4e+54 in total; families with more than one: 127 of 376
#> Two uniform draws differ on average in 84.6 of the 452 selected individuals (18.7%)

# A largest set with half women
rkn_select_global(g, "sex", target = c(F = 0.5), seed = 1)
#> <rkn_selection> Global selection (proportion of sex = F), seed 1
#> Call: rkn_select_global(g = g, covariate = "sex", target = c(F = 0.5), seed = 1)
#> Individuals: 751 | selected: 452 (60.2%)
#> Families: 376 | size min/mean/max: 1 / 2.00 / 31
#> Families with up to 5 members: 357 | with 6 or more: 19
#> Proportion of sex = F among the selected: 0.5000 (222 of 444 with a known value; target 0.5)
#> Selected with a missing value of sex: 8
#> Selections reaching this proportion: 1.5e+53, one drawn uniformly
```

## Features

- **Exact maximum sets**: every family is resolved by an exact
  algorithm, with a reduction that keeps intricate families tractable.
- **Selection uncertainty**: uniform draws among all the largest
  selections, their number, and the exact inclusion probability of every
  individual.
- **Preferences**: criteria on covariates, such as a minimum call rate
  or the age at onset of a disease, decide among the largest selections
  of each family.
- **Composition targets**: a target proportion or mean of a covariate
  over the whole selection, reached exactly by choosing all the families
  jointly, with the trade-off between selection size and composition.
- **Sample checks**: possible duplicates, contaminated samples and
  samples of different ancestry are flagged before selecting.
- **Diagnostics**: every decision is recorded and can be summarised and
  plotted.

## Documentation

- `vignette("RedKinNet")`: the whole workflow.
- `vignette("covariates")`: selecting with covariates.
- `vignette("uncertainty")`: quantifying the uncertainty of a selection.

## Citation

Run `citation("RedKinNet")` for the reference to cite.
