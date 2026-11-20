
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RedKinNet <img src="man/figures/logo.png" align="right" height="138" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/PaoloDalena/RedKinNet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaoloDalena/RedKinNet/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Many genetic analyses, such as association studies and principal component analysis, assume
that individuals are unrelated, so related samples are usually reduced to a set of unrelated
individuals. Existing tools either remove relatives greedily, which can discard more
individuals than necessary, or search exhaustively, which does not scale to the large families
of genetic isolates; and none of them says how much the result depends on which of the many
equally large sets is chosen.

RedKinNet selects unrelated individuals from pairwise kinship estimates, for instance those of
[KING](https://www.kingrelatedness.com/) or of
[PLINK 2](https://www.cog-genomics.org/plink/2.0/distance#make_king). Related pairs form a
network, which RedKinNet splits into families and solves exactly: the selection always keeps as
many individuals as possible, even in the large, intricate families of genetic isolates. Among the many
largest selections, RedKinNet can:

- draw one uniformly at random, and give the exact probability that each individual is selected;
- choose by preferences for individuals, such as a higher call rate;
- reach a target proportion or mean of a covariate over the whole selection, such as half women.

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
