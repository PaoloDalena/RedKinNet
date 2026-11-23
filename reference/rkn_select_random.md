# Select unrelated individuals at random

Selects a maximum set of unrelated individuals, a set with no related
pair of the largest possible size (see
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)),
and makes at random the choices that the network leaves open. The
network is partitioned into families, its connected components, and each
family is resolved on its own: since no relationship crosses two
families, the union of their selections has no related pair. By default,
each family is drawn uniformly at random among all its maximum sets by
an exact algorithm that also handles the intricate families of genetic
isolates (see the section Search over maximum sets of
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md));
alternatively, structural rules decide most individuals.

## Usage

``` r
rkn_select_random(g, uniform = TRUE, seed = NULL, diagnostics = TRUE)
```

## Arguments

- g:

  An igraph object. Vertex names must be stored in V(g)\$name. Vertex
  attributes cannot be named `degree`, `max_path_length` or `.name`:
  these names are reserved for the structural measures computed during
  the selection. The network is usually built with
  [`rkn_graph_from_kinship`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md),
  which also checks the samples; for a network built in another way, see
  [`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md).

- uniform:

  Logical. If `TRUE` (default), each family is drawn uniformly at random
  among all its maximum sets, so that every maximum set of the network
  has the same probability, and the selection probabilities are exact:
  with diagnostics, `node_diagnostics$inclusion_probability` gives the
  probability that each individual is selected. Families too intricate
  for this exact computation are selected as with `uniform = FALSE` and
  reported in the column `uniform` of `family_summary` (see
  [`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)).
  If `FALSE`, structural rules decide most individuals and the remaining
  choices are made at random; the result is still a maximum set, but the
  maximum sets are not equally likely, and the rules never choose some
  of them: in a chain of six relatives, two of its four maximum sets.
  The rules record why each individual was kept or discarded, and are
  somewhat faster in cohorts whose families are small.

- seed:

  Optional integer. If provided, the selection is reproducible: a seed
  is drawn from it for each family, so that the selection within a
  family does not depend on the random numbers used by the others, and
  different values give independent selections. If `NULL` (default), the
  current random number stream is used.

- diagnostics:

  Logical. If `TRUE` (default), the returned object also records how
  each individual and each family was resolved (`node_diagnostics` and
  `family_summary`), which
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  `plot(type = "diagnostics")` use. Diagnostics increase the run time
  and the size of the returned object: computing the inclusion
  probabilities makes a selection two to three times slower in family
  studies, though hardly slower in population cohorts. On very large
  networks, or when many selections are drawn, `diagnostics = FALSE`
  returns the selection alone.

## Value

An object of class
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md),
with `print`, `summary` and `plot` methods; `$kept` holds the
identifiers of the selected individuals.

## Computation time

On one core of a laptop (Apple M2), the uniform selection with
diagnostics took about 0.3 seconds for a population cohort of 1,000
samples, 3 seconds for 10,000 and 30 seconds for 100,000, and the
selection by rules about two thirds of these times. The time depends on
the largest and most intricate families more than on the number of
samples: a family study of 10,000 samples took 40 seconds (16 without
diagnostics), and a family study of 100,000 samples about 12 minutes on
one core of a compute node, where a population cohort of the same size
took less than 2 minutes.

## See also

[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md),
[`summary.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/summary.rkn_selection.md),
[`plot.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md),
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md),
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md).

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
# Uniformly among all the maximum sets, with the probability that each individual is
# selected
sel <- rkn_select_random(g, seed = 1)
sel
#> <rkn_selection> Uniform random selection, seed 1
#> Call: rkn_select_random(g = g, seed = 1)
#> Individuals: 15 | selected: 8 (53.3%)
#> Families: 5 | size min/mean/max: 1 / 3.00 / 6
#> Families with up to 5 members: 4 | with 6 or more: 1
#> Families sampled uniformly among their maximum sets: 5 of 5
#> Maximum sets: 6 in total; families with more than one: 2 of 5
#> Two uniform draws differ on average in 1.2 of the 8 selected individuals (14.6%)
sel$kept
#> [1] "P1" "P2" "S2" "T3" "D1" "D2" "D3" "U1"
sel$node_diagnostics[, c("id", "inclusion_probability")]
#>    id inclusion_probability
#> 1  P1             1.0000000
#> 2  C1             0.0000000
#> 3  P2             1.0000000
#> 4  S1             0.5000000
#> 5  S2             0.5000000
#> 6  T1             0.3333333
#> 7  T2             0.3333333
#> 8  T3             0.3333333
#> 9  D1             1.0000000
#> 10 K1             0.0000000
#> 11 D2             1.0000000
#> 12 K2             0.0000000
#> 13 K3             0.0000000
#> 14 D3             1.0000000
#> 15 U1             1.0000000

# The selection alone, without diagnostics
rkn_select_random(g, seed = 1, diagnostics = FALSE)$kept
#> [1] "P1" "P2" "S2" "T3" "D1" "D2" "D3" "U1"

# By structural rules, which record why each individual was kept or discarded
rules <- rkn_select_random(g, uniform = FALSE, seed = 1)
summary(rules)
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
```
