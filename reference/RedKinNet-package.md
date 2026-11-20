# RedKinNet: selection of unrelated individuals from kinship networks

Genetic studies often need a set of individuals with no close relatives
among them. RedKinNet joins every pair of related individuals in a
network and selects a maximum set of unrelated individuals: a set with
no related pair of the largest possible size (see
[`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)).
The families, the connected components of the network, are resolved by
structural rules and by an exact search that also handles the large and
intricate families of genetic isolates.

A typical analysis follows these steps:

1.  [`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md)
    flags possible duplicates and contaminated samples in the kinship
    estimates;

2.  [`rkn_graph_from_kinship`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md)
    builds the network at a kinship threshold, optionally leaving
    samples out and adding covariates;

3.  a selector chooses the individuals:
    [`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
    at random, optionally uniformly among all the maximum sets with
    exact inclusion probabilities;
    [`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
    by preferences for individuals, such as women; or
    [`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
    towards a target composition of the whole selection, such as half
    women or a mean age, after
    [`rkn_tradeoff`](https://paolodalena.github.io/RedKinNet/reference/rkn_tradeoff.md)
    has shown what the target costs in selection size;

4.  the returned
    [`rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_selection.md)
    object is inspected with
    [`print()`](https://rdrr.io/r/base/print.html),
    [`summary()`](https://rdrr.io/r/base/summary.html) and
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (see
    [`rkn_plot_pre`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md)).

## See also

Useful links:

- <https://github.com/PaoloDalena/RedKinNet>

- <https://paolodalena.github.io/RedKinNet/>

- Report bugs at <https://github.com/PaoloDalena/RedKinNet/issues>

## Author

**Maintainer**: Paolo Dalena <paolodalena97@gmail.com>
([ORCID](https://orcid.org/0000-0003-2920-9572)) \[copyright holder\]
