# Prefer individuals above or below a cutoff

Scores a candidate set of individuals by the proportion of its members
whose `covariate` value lies on the preferred side of `cutoff`.

Criteria are used by
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md):
within each family, where several sets of unrelated individuals remain
to choose from, each candidate set is scored and the best one is
selected.

## Usage

``` r
rkn_criterion_cutoff(covariate, cutoff, prefer = c("above", "below"))
```

## Arguments

- covariate:

  Character. Name of a numeric vertex attribute of the network, such as
  `"age"` for `V(g)$age`.

- cutoff:

  Numeric. Threshold value.

- prefer:

  Character. Either `"above"` (prefer `covariate > cutoff`) or `"below"`
  (prefer `covariate < cutoff`).

## Value

An object of class `c("rkn_criterion_cutoff", "rkn_criterion")`.

## See also

[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md),
[`rkn_criterion_level`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_level.md);
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
for a target composition of the whole selection.

## Examples

``` r
# Prefer individuals younger than 50
rkn_criterion_cutoff("age", cutoff = 50, prefer = "below")
#> <rkn_criterion: cutoff>
#>   covariate: age
#>   cutoff: 50
#>   prefer: below
```
