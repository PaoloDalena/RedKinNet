# RedKinNet 0.1.0

First release.

* `rkn_graph_from_kinship()` builds the network of related samples from a kinship matrix or
  a table of pairs, with thresholds given as numbers or as KING degrees of relationship.
* `rkn_check_samples()` flags possible duplicates, contaminated samples and samples of
  different ancestry.
* `rkn_select_random()` selects a largest set of unrelated individuals by structural rules
  or, with `uniform = TRUE`, uniformly at random among all of them, with exact inclusion
  probabilities.
* `rkn_select_criteria()` chooses among the largest sets by preferences for individuals
  (`rkn_criterion_cutoff()`, `rkn_criterion_level()`), within the choices left by the rules
  or, with `optimal = TRUE`, over all the largest sets of each family.
* `rkn_select_global()` reaches a target proportion or mean of a covariate over the whole
  selection exactly, optionally giving up individuals (`max_loss`); `rkn_tradeoff()` shows
  the composition reachable for each number of individuals given up.
* Selections are objects of class `rkn_selection` with `print()`, `summary()` and `plot()`
  methods; `rkn_plot_pre()` draws a network before the selection.
* `rkn_example` is a synthetic cohort for the examples and the three vignettes.
