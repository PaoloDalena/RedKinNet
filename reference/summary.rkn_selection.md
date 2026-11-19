# Summarise an rkn_selection object

Computes summary tables of a selection of unrelated individuals: the
selection by family structure and, when diagnostics were computed, how
each individual was decided. The family structures and the types of
decision are described in the sections Family structures and Types of
decision of
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md).

## Usage

``` r
# S3 method for class 'rkn_selection'
summary(object, ...)
```

## Arguments

- object:

  An `rkn_selection` object, as returned by
  [`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md),
  [`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
  or
  [`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md).

- ...:

  Not used.

## Value

An object of class `summary.rkn_selection`: a list with elements

- `header`:

  a one-line description of how the selection was made;

- `overall`:

  a one-row data frame with the numbers of individuals, of selected
  individuals and of families;

- `by_structure`:

  a data frame with one row per family structure present, giving the
  family size (a range for families of 6 or more members), the
  structure, and the numbers of families, of individuals and of selected
  individuals;

- `decisions`:

  a contingency table of individuals by type of decision and outcome
  (selected or discarded), or `NULL` if diagnostics were not computed
  and for global selections, which have no types of decision;

- `method`:

  the `method` of the selection.

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
summary(sel)
#> <rkn_selection summary> Random selection, seed 1
#> Individuals: 15 | selected: 8 (53.3%) | families: 5
#> 
#> Selection by family structure:
#>  size                               structure families individuals selected
#>     1                       single individual        1           1        1
#>     2                            related pair        1           2        1
#>     3 chain of three (e.g. parents and child)        1           3        2
#>     3                  three mutually related        1           3        1
#>     6             family of 6 or more members        1           6        3
#> 
#> Individuals by type of decision and outcome:
#>                                  outcome
#> decision                          selected discarded
#>   kept without a choice                  6         0
#>   discarded without a choice             0         4
#>   decided among single candidates        2         3
#>   decided among alternative sets         0         0

# The tables as data
s <- summary(sel)
s$by_structure
#>   family_size                               structure n_families n_individuals
#> 1           1                       single individual          1             1
#> 2           2                            related pair          1             2
#> 3           3 chain of three (e.g. parents and child)          1             3
#> 4           3                  three mutually related          1             3
#> 5           6             family of 6 or more members          1             6
#>   n_selected
#> 1          1
#> 2          1
#> 3          2
#> 4          1
#> 5          3
s$decisions
#>                                  outcome
#> decision                          selected discarded
#>   kept without a choice                  6         0
#>   discarded without a choice             0         4
#>   decided among single candidates        2         3
#>   decided among alternative sets         0         0
```
