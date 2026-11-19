# Package index

## Kinship network

Build the network of related samples and check the samples.

- [`rkn_graph_from_kinship()`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md)
  : Build a relatedness network from kinship estimates
- [`rkn_check_samples()`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md)
  : Check samples for duplicates, possible contamination and different
  ancestry
- [`print(`*`<rkn_sample_check>`*`)`](https://paolodalena.github.io/RedKinNet/reference/print.rkn_sample_check.md)
  : Print an rkn_sample_check object

## Selection of unrelated individuals

Select a largest set of unrelated individuals at random, by preferences
for individuals, or with a target composition.

- [`rkn_select_random()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
  : Select unrelated individuals at random
- [`rkn_select_criteria()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
  : Select unrelated individuals by criteria
- [`rkn_criterion_cutoff()`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_cutoff.md)
  : Prefer individuals above or below a cutoff
- [`rkn_criterion_level()`](https://paolodalena.github.io/RedKinNet/reference/rkn_criterion_level.md)
  : Prefer individuals with a given value of a covariate
- [`rkn_select_global()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
  : Select unrelated individuals with a target composition
- [`rkn_tradeoff()`](https://paolodalena.github.io/RedKinNet/reference/rkn_tradeoff.md)
  [`print(`*`<rkn_tradeoff>`*`)`](https://paolodalena.github.io/RedKinNet/reference/rkn_tradeoff.md)
  [`plot(`*`<rkn_tradeoff>`*`)`](https://paolodalena.github.io/RedKinNet/reference/rkn_tradeoff.md)
  : Compute the trade-off between selection size and composition

## Selections

Results of the selectors, their summaries and plots.

- [`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)
  : Selections of unrelated individuals
- [`print(`*`<rkn_selection>`*`)`](https://paolodalena.github.io/RedKinNet/reference/print.rkn_selection.md)
  : Print an rkn_selection object
- [`summary(`*`<rkn_selection>`*`)`](https://paolodalena.github.io/RedKinNet/reference/summary.rkn_selection.md)
  : Summarise an rkn_selection object
- [`print(`*`<summary.rkn_selection>`*`)`](https://paolodalena.github.io/RedKinNet/reference/print.summary.rkn_selection.md)
  : Print a summary of an rkn_selection object
- [`plot(`*`<rkn_selection>`*`)`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md)
  [`rkn_plot_pre()`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md)
  : Plot a relatedness network

## Data

- [`rkn_example`](https://paolodalena.github.io/RedKinNet/reference/rkn_example.md)
  : Example cohort with simulated kinship estimates
- [`rules_db`](https://paolodalena.github.io/RedKinNet/reference/rules_db.md)
  : Structural rules for families of up to 5 members

## Package

- [`RedKinNet`](https://paolodalena.github.io/RedKinNet/reference/RedKinNet-package.md)
  [`RedKinNet-package`](https://paolodalena.github.io/RedKinNet/reference/RedKinNet-package.md)
  : RedKinNet: selection of unrelated individuals from kinship networks
