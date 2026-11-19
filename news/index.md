# Changelog

## RedKinNet 0.1.0

First release.

- [`rkn_graph_from_kinship()`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md)
  builds the network of related samples from a kinship matrix or a table
  of pairs, with thresholds given as numbers or as KING degrees of
  relationship.
- [`rkn_check_samples()`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md)
  flags possible duplicates, contaminated samples and samples of
  different ancestry.
- [`rkn_select_random()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
  selects a largest set of unrelated individuals by structural rules or,
  with `uniform = TRUE`, uniformly at random among all of them, with
  exact inclusion probabilities.
- [`rkn_select_criteria()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
  chooses among the largest sets by preferences for individuals
  ([`rkn_criterion_cutoff()`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_cutoff.md),
  [`rkn_criterion_level()`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_level.md)),
  within the choices left by the rules or, with `optimal = TRUE`, over
  all the largest sets of each family.
- [`rkn_select_global()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
  reaches a target proportion or mean of a covariate over the whole
  selection exactly, optionally giving up individuals (`max_loss`);
  [`rkn_tradeoff()`](https://paolodalena.github.io/RedKinNet/reference/rkn_tradeoff.md)
  shows the composition reachable for each number of individuals given
  up.
- Selections are objects of class `rkn_selection` with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods;
  [`rkn_plot_pre()`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md)
  draws a network before the selection.
- `rkn_example` is a synthetic cohort for the examples and the three
  vignettes.
