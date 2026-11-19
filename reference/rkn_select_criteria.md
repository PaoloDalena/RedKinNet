# Select unrelated individuals by criteria

Selects a maximum set of unrelated individuals, as
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
does, but makes the choices that the network leaves open by criteria on
the covariates instead of at random. Each family is resolved on its own:
families of up to 5 members through the structural rules in
[`rules_db`](https://paolodalena.github.io/RedKinNet/reference/rules_db.md),
larger families through structural reductions followed by a search of
what remains; the returned object is as for
[`rkn_select_random()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md).
How the criteria are evaluated in large families, exactly or on a sample
of maximum sets, is described in the section Search over maximum sets of
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md).
With `optimal = TRUE`, the structural rules are skipped and the criteria
choose among all the maximum sets of each family.

Covariates referenced by the criteria are read from the vertex
attributes of `g`. The criteria express preferences for individuals,
such as women or people younger than 50; for a target composition of the
whole selection, such as half women or a mean age, see
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md).

## Usage

``` r
rkn_select_criteria(
  g,
  criterion_primary,
  criterion_secondary = NULL,
  optimal = FALSE,
  seed = NULL,
  diagnostics = TRUE
)
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

- criterion_primary:

  A criterion, created by
  [`rkn_criterion_level`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_level.md)
  or
  [`rkn_criterion_cutoff`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_cutoff.md),
  used to score candidate resolutions within each family.

- criterion_secondary:

  Optional `rkn_criterion` object used to break ties in
  `criterion_primary`.

- optimal:

  Logical. If `FALSE` (default), structural rules decide part of each
  family before the criteria are applied, and some of them fix
  individuals regardless of the criteria (for instance, they keep the
  individual at the end of a chain of relatives); the result is the best
  maximum set among those compatible with the rules. If `TRUE`, the
  criteria are evaluated over all the maximum sets of each family,
  exactly, or as described in the section Search over maximum sets of
  [`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)
  when the covariate has missing values. Families too intricate for this
  computation are selected as with `optimal = FALSE` and reported in the
  column `optimal` of `family_summary`.

- seed:

  Optional integer, as in
  [`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md):
  a seed is drawn from it for each family. Random numbers are used only
  to break exact ties in the criterion scores.

- diagnostics:

  Logical. If `TRUE` (default), the returned object also records how
  each individual and each family was resolved (`node_diagnostics` and
  `family_summary`), which
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  `plot(type = "diagnostics")` use. Diagnostics increase the run time
  and the size of the returned object; on very large networks,
  `diagnostics = FALSE` returns the selection alone.

## Value

An object of class
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md),
as for
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md);
its `criteria` element holds the criteria.

## See also

[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md),
[`summary.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/summary.rkn_selection.md),
[`plot.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md),
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md),
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md),
[`rkn_criterion_level`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_level.md),
[`rkn_criterion_cutoff`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_cutoff.md),
[`rules_db`](https://paolodalena.github.io/RedKinNet/reference/rules_db.md).

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
# Where the structure of a family leaves a choice, prefer women ...
rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1)$kept
#> [1] "P1" "P2" "S2" "T2" "D3" "D1" "D2" "U1"

# ... evaluating the criterion over all the maximum sets of each family ...
rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1, optimal = TRUE)$kept
#> [1] "P1" "P2" "S2" "T2" "D1" "D2" "D3" "U1"

# ... or prefer individuals younger than 40
rkn_select_criteria(g, rkn_criterion_cutoff("age", 40, prefer = "below"), seed = 1)$kept
#> [1] "P1" "P2" "S2" "T2" "D3" "D1" "D2" "U1"

# A secondary criterion breaks the ties left by the primary one
sel <- rkn_select_criteria(g, rkn_criterion_level("sex", "F"),
                           rkn_criterion_cutoff("age", 40, prefer = "below"), seed = 1)
sel
#> <rkn_selection> Criteria-based selection (level, then cutoff), seed 1
#> Call: rkn_select_criteria(g = g, criterion_primary = rkn_criterion_level("sex", "F"), criterion_secondary = rkn_criterion_cutoff("age", 40, prefer = "below"), seed = 1)
#> Individuals: 15 | selected: 8 (53.3%)
#> Families: 5 | size min/mean/max: 1 / 3.00 / 6
#> Families with up to 5 members: 4 | with 6 or more: 1
#> Families resolved by a search over maximum sets: 1
```
