# Selections of unrelated individuals

[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md),
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
and
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
return an object of class `rkn_selection`: a list holding the selected
individuals, the network they were selected from and, when computed,
diagnostics on how each individual and each family was resolved. The
object has
[`print`](https://paolodalena.github.io/RedKinNet/reference/print.rkn_selection.md),
[`summary`](https://paolodalena.github.io/RedKinNet/reference/summary.rkn_selection.md)
and
[`plot`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md)
methods.

Throughout the package, a *maximum set* of unrelated individuals is a
set with no related pair of the largest possible size: in graph terms, a
maximum independent set of the relatedness network. A network usually
has many maximum sets, which differ in which individuals they contain; a
*family* is a connected component of the network (see the section Family
structures).

## Components

- `kept`:

  character vector, the identifiers of the selected individuals;

- `discarded`:

  character vector, the identifiers of the other individuals;

- `family_sizes`:

  integer vector, the size of each family (connected component of the
  network);

- `method`:

  `"random"`, `"criteria"` or `"global"`;

- `uniform`:

  logical, whether the selection was drawn uniformly among the maximum
  sets (`rkn_select_random(uniform = TRUE)`);

- `optimal`:

  logical, whether the criteria were evaluated over all the maximum sets
  of each family (`rkn_select_criteria(optimal = TRUE)`);

- `criteria`:

  for criteria-based selections, a list with the `primary` and, if
  given, the `secondary` criterion; `NULL` otherwise;

- `seed`:

  the seed of the selection, or `NULL`;

- `diagnostics`:

  logical, whether diagnostics were computed;

- `node_diagnostics`:

  a data frame with one row per individual, or `NULL` without
  diagnostics. Columns: `id`; `pre_class`, how the individual was
  decided (`"kept_apriori"`, `"discarded_apriori"`, `"candidate_simple"`
  or `"candidate_complex"`, described with their labels in the section
  Types of decision below); `post_class`, `"kept"` or `"discarded"`;
  `is_selected`, the outcome as a logical; `pre_reason`, an internal
  label of the rule or procedure that decided the individual;
  `family_id` and `family_size`; in uniform selections,
  `inclusion_probability`, the exact probability that the individual is
  selected, i.e. the proportion of the maximum sets of its family that
  contain it (`NA` in the families too intricate for uniform sampling).
  From these probabilities
  [`print()`](https://rdrr.io/r/base/print.html) reports the expected
  number of selected individuals that differ between two independent
  uniform draws, the sum of \\p(1 - p)\\ over the individuals: the part
  of the selection left to chance;

- `family_summary`:

  a data frame with one row per family, or `NULL` without diagnostics.
  Columns: `family_id` and `family_size`; `selector`, `"upto5"` for
  families resolved by the structural rules in
  [`rules_db`](https://paolodalena.github.io/RedKinNet/reference/rules_db.md),
  `"6plus"` for larger families, `"uniform"` for families drawn
  uniformly and `"optimal"` for families where the criteria were
  evaluated over all maximum sets; `structure_name`, an internal label
  of the recognised structure; `mis_method`, how the individuals left
  undecided by the structural rules were resolved (see the section
  Search over maximum sets): `"enumeration"`, `"exact"`,
  `"exact_reduced"` or `"exact_sample"`, or `NA` when no search was
  needed because the rules decided every individual or chose among the
  candidates they identified, as in a related pair; in uniform
  selections and in selections over all maximum sets, every family is
  searched; `uniform`, whether the family was drawn uniformly among its
  maximum sets; `optimal`, whether the criteria were evaluated over all
  its maximum sets; `n_max_sets`, in the families drawn uniformly or
  evaluated over all maximum sets, the number of maximum sets of the
  family (`NA` otherwise, and after `"exact_reduced"`; `Inf` beyond
  about \\1.8 \times 10^{308}\\, the largest number R can store); the
  numbers of kept and discarded individuals, overall (`n_kept`,
  `n_discarded`) and by type of decision (`n_kept_apriori`,
  `n_discarded_apriori`, `n_candidate_simple`, `n_candidate_complex`).
  Since the choices in different families are independent, the maximum
  sets of the network are all the combinations of one maximum set per
  family: their number, reported by
  [`print()`](https://rdrr.io/r/base/print.html), is the product of
  `n_max_sets` over the families and grows exponentially with the number
  of families that have a choice (each related pair doubles it), so
  [`print()`](https://rdrr.io/r/base/print.html) also reports how many
  families have more than one maximum set. In global selections the
  columns are `family_id`, `family_size`, `selector` (`"global"`),
  `exact` (`FALSE` for the families beyond the computational budget, see
  [`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)),
  `n_kept`, `n_discarded`, `n_category` (selected members in the
  category) or, for a numeric covariate, `mean` (mean of the covariate
  among the selected members with a known value, `NA` if none), and
  `n_missing` (selected members with a missing value), and the types of
  decision of `node_diagnostics` are not available;

- `graph`:

  the igraph network the selection was made on;

- `call`:

  the call that produced the selection;

- `global`:

  for global selections, a list describing the composition reached:
  `type`, `"proportion"` for a categorical covariate or `"mean"` for a
  numeric one; `covariate` and `target`; `cohort`, whether the target is
  the proportion or mean of the network; `n_known` and `n_missing`, the
  selected individuals with a known and a missing value; `size`,
  `max_size` and `max_loss`; `log_count`, the natural logarithm of the
  number of selections reaching the same composition and size (`NA` when
  `exact` is `FALSE`); `exact`, whether every family was computed
  exactly; `n_inexact`, the number of families that were not. For a
  categorical covariate, also `category`, `proportion` (the proportion
  reached among the selected individuals with a known value) and
  `n_category`; for a numeric covariate, `mean` (the mean of the
  recorded values among the selected individuals with a known value),
  `resolution`, `rounded` (whether some value was rounded for the
  search) and `coarsened` (whether the default resolution was coarsened
  to fit the memory and time budgets). `NULL` otherwise.

## Family structures

A family is a connected component of the relatedness network: a group of
individuals linked, directly or through other members, by relatedness
above the threshold used to build the network.
[`summary.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/summary.rkn_selection.md)
groups the families into seven structures:

- single individual:

  an individual related to no one else in the cohort;

- related pair:

  two related individuals;

- chain of three:

  three individuals of whom two are not related to each other but both
  are related to the third, typically two parents and their child;

- three mutually related:

  three individuals all related to each other, for instance three
  siblings or a parent and two children;

- family of 4, of 5, and of 6 or more members:

  larger families. Families of up to 5 members are resolved by the
  structural rules in
  [`rules_db`](https://paolodalena.github.io/RedKinNet/reference/rules_db.md),
  larger families by reducing the family step by step with structural
  rules and resolving what remains.

## Search over maximum sets

The structural rules decide many individuals directly. The individuals
they leave undecided are resolved by a search over the maximum sets of
unrelated individuals, reported in the column `mis_method` of
`family_summary`:

- `"enumeration"`:

  in families of up to 5 members, all the maximum sets are listed and
  one of them is chosen;

- `"exact"`:

  in families of 6 or more members, an exact algorithm (dynamic
  programming along an elimination ordering of the family) finds the
  maximum size and draws one of the maximum sets uniformly at random.
  Its cost grows with how intricate the relationships are, not with the
  size of the family, so it also handles the large families of genetic
  isolates;

- `"exact_reduced"`:

  the same after removing the individuals who cannot enlarge the
  selection (an individual is removed when one of their relatives has no
  relative outside their own). This is used only for families too
  intricate for the exact algorithm: the set chosen is still of maximum
  size, but not drawn uniformly. In criteria-based selections an
  individual is removed only if that relative meets the criteria at
  least as well, so the criteria are still optimised exactly, or, with
  missing covariate values, evaluated on the maximum sets that remain.
  If the family is still too intricate, individuals are removed
  regardless of the criteria, as in selections without criteria: the set
  is still of maximum size, and the criteria are optimised over the
  maximum sets that remain.

- `"exact_sample"`:

  criteria-based selections only: see below.

In criteria-based selections
([`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)),
the criteria are optimised exactly by the same algorithm when no
undecided member of the family lacks the covariate, because among sets
of the same size they rank sets by the number of members meeting the
condition. When some undecided member lacks it, the score is a
proportion among the members with a known value, which the algorithm
cannot maximise directly: the criteria are then evaluated on all the
maximum sets when there are at most 10,000; above that number, as in
some large families of genetic isolates, they are evaluated on a uniform
random sample of 10,000 maximum sets, and the best of the sample is kept
(`"exact_sample"`). The number of such families is reported by
[`print()`](https://rdrr.io/r/base/print.html). In uniform selections
and in selections over all maximum sets, no structural rule is applied,
and every family is searched in this way. If a family remains too
intricate even after this reduction, the selection stops with an error
that names one of its members and suggests excluding the samples flagged
by
[`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md)
or using a higher threshold, which joins closer relatives only.

## Types of decision

Within each family, every individual is decided in one of four ways,
reported by
[`summary.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/summary.rkn_selection.md)
and shown by the `"diagnostics"` plot (see
[`plot.rkn_selection`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md)):

- kept without a choice (`"kept_apriori"`):

  the individual is kept with no alternative to choose from: a rule
  based on the structure of the family keeps them (for instance the
  leaves of a star-shaped family), or they belong to every maximum set
  of the individuals left undecided by the rules;

- discarded without a choice (`"discarded_apriori"`):

  likewise, the individual is discarded by a structural rule (for
  instance the centre of a star) or belongs to no maximum set of the
  individuals left undecided;

- decided among single candidates (`"candidate_simple"`):

  the remaining choice is among interchangeable candidates of a single
  group, such as the two members of a related pair, and one of them is
  selected;

- decided among alternative sets (`"candidate_complex"`):

  several maximum sets remain that differ in more than one individual,
  and one of them is selected. These are the sets left by the rules, not
  necessarily all the maximum sets of the family.

The choices among candidates are made at random by
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
and by the criteria in
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md).

"Without a choice" describes the procedure, not the family: to simplify
the search, the structural rules fix some individuals that other maximum
sets of the family leave out (for instance, in a chain of relatives they
keep the individual at the end), so an individual kept without a choice
is not necessarily in every maximum set of the family. The exact classes
over the whole family are given by uniform selections.

In uniform selections (`rkn_select_random(uniform = TRUE)`) no
structural rule is applied, and the same codes have an exact meaning,
computed from the inclusion probabilities: `"kept_apriori"`, in every
maximum set (probability 1); `"discarded_apriori"`, in no maximum set
(probability 0); `"candidate_complex"`, chosen at random among the
maximum sets (probability strictly between 0 and 1).
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) label them
accordingly.

In criteria-based selections over all maximum sets
(`rkn_select_criteria(optimal = TRUE)`) no structural rule is applied
either: `"kept_apriori"` and `"discarded_apriori"` again mean in every
and in no maximum set, and `"candidate_complex"` chosen by the criteria
among maximum sets that differ, which
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) label
accordingly.

In both kinds of selection, the families too intricate for the exact
computation, whose number [`print()`](https://rdrr.io/r/base/print.html)
reports, are selected by the structural rules, and their codes keep the
meaning given above for selections with rules. In the families reduced
by dominance (`mis_method` `"exact_reduced"`), the individuals removed
by the reduction count as in no maximum set, although some maximum sets
may contain them.

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
names(sel)
#>  [1] "kept"             "discarded"        "family_sizes"     "method"          
#>  [5] "uniform"          "optimal"          "criteria"         "seed"            
#>  [9] "diagnostics"      "node_diagnostics" "family_summary"   "graph"           
#> [13] "call"             "global"          
sel$kept
#> [1] "P1" "P2" "S2" "T3" "D1" "D2" "D3" "U1"
head(sel$node_diagnostics)
#>   id         pre_class post_class is_selected     pre_reason
#> 1 P1      kept_apriori       kept        TRUE  uniform_every
#> 2 C1 discarded_apriori  discarded       FALSE   uniform_none
#> 3 P2      kept_apriori       kept        TRUE  uniform_every
#> 4 S1 candidate_complex  discarded       FALSE uniform_chosen
#> 5 S2 candidate_complex       kept        TRUE uniform_chosen
#> 6 T1 candidate_complex  discarded       FALSE uniform_chosen
#>   inclusion_probability family_id family_size
#> 1             1.0000000         1           3
#> 2             0.0000000         1           3
#> 3             1.0000000         1           3
#> 4             0.5000000         2           2
#> 5             0.5000000         2           2
#> 6             0.3333333         3           3
sel$family_summary[, c("family_id", "family_size", "mis_method", "n_kept", "n_discarded")]
#>   family_id family_size mis_method n_kept n_discarded
#> 1         1           3      exact      2           1
#> 2         2           2      exact      1           1
#> 3         3           3      exact      1           2
#> 4         4           6      exact      3           3
#> 5         5           1      exact      1           0
```
