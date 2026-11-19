# Print an rkn_selection object

Prints a brief summary of a selection of unrelated individuals: how the
selection was made, how many individuals were selected, and how the
network decomposes into families.

## Usage

``` r
# S3 method for class 'rkn_selection'
print(x, ...)
```

## Arguments

- x:

  An `rkn_selection` object, as returned by
  [`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md),
  [`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
  or
  [`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md).

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)

## Examples

``` r
# A small relatedness network: two parents and their child (P1, P2, C1), a pair of
# siblings (S1, S2), three mutually related individuals (T1, T2, T3), a larger family
# (parents D1 and D2, their children K1 and K2, and K3, born to D2 and D3, a
# half-sibling of K1 and K2) and an unrelated individual (U1).
g <- igraph::graph_from_literal(P1 - C1 - P2, S1 - S2, T1 - T2 - T3 - T1,
                                D1 - K1 - D2, D1 - K2 - D2, K1 - K2,
                                D2 - K3 - D3, K1 - K3, K2 - K3, U1)
sex <- c(P1 = "F", C1 = "M", P2 = "M", S1 = "M", S2 = "F", T1 = "M", T2 = "F", T3 = "M",
         D1 = "M", K1 = "F", D2 = "F", K2 = "M", K3 = "F", D3 = "M", U1 = "F")
age <- c(P1 = 52, C1 = 21, P2 = 55, S1 = 34, S2 = 31, T1 = 60, T2 = 33, T3 = 29,
         D1 = 63, K1 = 35, D2 = 61, K2 = 38, K3 = 30, D3 = 58, U1 = 44)
igraph::V(g)$sex <- unname(sex[igraph::V(g)$name])
igraph::V(g)$age <- unname(age[igraph::V(g)$name])
sel <- rkn_select_random(g, seed = 1)
sel # the same as print(sel)
#> <rkn_selection> Random selection, seed 1
#> Call: rkn_select_random(g = g, seed = 1)
#> Individuals: 15 | selected: 8 (53.3%)
#> Families: 5 | size min/mean/max: 1 / 3.00 / 6
#> Families with up to 5 members: 4 | with 6 or more: 1
#> Families resolved by a search over maximum sets: 1

print(rkn_select_random(g, seed = 1, diagnostics = FALSE))
#> <rkn_selection> Random selection, seed 1
#> Call: rkn_select_random(g = g, seed = 1, diagnostics = FALSE)
#> Individuals: 15 | selected: 8 (53.3%)
#> Families: 5 | size min/mean/max: 1 / 3.00 / 6
#> Families with up to 5 members: 4 | with 6 or more: 1
#> Node-level diagnostics were not computed (diagnostics = FALSE); rerun with diagnostics = TRUE for a richer output.
```
