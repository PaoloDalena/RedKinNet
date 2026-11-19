# Prefer individuals with a given value of a covariate

Scores a candidate set of individuals by the proportion of its members
whose `covariate` value equals `level`. Intended for dichotomous
covariates, but works for any covariate type supporting `==`.

Criteria are used by
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md):
within each family, where several sets of unrelated individuals remain
to choose from, each candidate set is scored and the best one is
selected.

## Usage

``` r
rkn_criterion_level(covariate, level)
```

## Arguments

- covariate:

  Character. Name of a vertex attribute of the network, such as `"sex"`
  for `V(g)$sex`.

- level:

  The preferred value of `covariate` (matched with `==`).

## Value

An object of class `c("rkn_criterion_level", "rkn_criterion")`.

## See also

[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md),
[`rkn_criterion_cutoff`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_cutoff.md);
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
for a target composition of the whole selection.

## Examples

``` r
# Prefer women
rkn_criterion_level("sex", "F")
#> <rkn_criterion: level>
#>   covariate: sex
#>   level: F
```
