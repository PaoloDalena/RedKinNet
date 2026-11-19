# Build a relatedness network from kinship estimates

Builds the relatedness network used by
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md),
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md),
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
and
[`rkn_tradeoff`](https://paolodalena.github.io/RedKinNet/reference/rkn_tradeoff.md):
one vertex per sample, and an edge between two samples whose kinship
coefficient is at least `threshold`. The families of the selectors are
the connected components of this network.

## Usage

``` r
rkn_graph_from_kinship(
  x,
  threshold,
  samples = NULL,
  exclude = NULL,
  covariates = NULL,
  covariates_id = "id",
  check = TRUE
)
```

## Arguments

- x:

  Kinship estimates, either a square symmetric matrix with the sample
  identifiers as row and column names (the diagonal is ignored), or a
  data frame with one row per pair of samples, such as the `.kin0`
  output of KING or the output of `plink2 --make-king-table`. In a data
  frame, the two identifier columns are recognised by the names `ID1`
  and `ID2` (KING), `IID1` and `IID2` (PLINK 2) or `id1` and `id2`, and
  the kinship column by the name `kinship`, regardless of case.

- threshold:

  The kinship coefficient from which two samples are considered related:
  a number, or one of `"duplicate"`, `"first"`, `"second"` and
  `"third"`, which select the cutoff of KING for relatives of that
  degree or closer (see Details). There is no default: the choice
  depends on the study.

- samples:

  Optional character vector with the identifiers of all genotyped
  samples, for instance the IID column of the PLINK `.fam` or `.psam`
  file. Pair tables usually list only related pairs, so samples without
  relatives appear in no row; `samples` adds them to the network as
  unrelated individuals. Ignored for a matrix.

- exclude:

  Optional character vector of sample identifiers to leave out of the
  network, for instance samples flagged by a quality check.

- covariates:

  Optional data frame of sample covariates, with the sample identifiers
  in the column named by `covariates_id`. Its other columns become
  vertex attributes of the network, which
  [`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
  and
  [`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md)
  read; samples missing from `covariates` get `NA`. Covariates cannot be
  named `name`, which holds the sample identifiers, nor `degree`,
  `max_path_length` or `.name`, which are reserved for the structural
  measures computed during the selection.

- covariates_id:

  Name of the identifier column of `covariates`.

- check:

  Logical. If `TRUE` (default), the samples left after `exclude` are
  checked with
  [`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md),
  with the same threshold and default settings, and a single warning
  summarises the possible duplicates and contamination it finds. Nothing
  is removed. Samples of possibly different ancestry, which do not
  affect the selection, are reported only by
  [`rkn_check_samples()`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md).
  `check = FALSE` skips the check, for instance when rebuilding a
  network whose samples were already checked.

## Value

An undirected igraph object with one vertex per sample: the vertex
attribute `name` holds the sample identifiers and the other vertex
attributes the covariates, the edge attribute `kinship` holds the
kinship coefficient, and the graph attribute `kinship_threshold` holds
the cutoff used.

## Details

The threshold labels select the boundaries between the degrees of
relationship inferred by KING (Manichaikul et al., 2010):

|  |  |  |
|----|----|----|
| **label** | **cutoff** | **relatives joined, for instance** |
| `"duplicate"` | 0.354 | the same individual genotyped twice, monozygotic twins |
| `"first"` | 0.177 | parent and child, full siblings, and the above |
| `"second"` | 0.0884 | grandparent and grandchild, aunt or uncle and nephew or niece, half-siblings, and the above |
| `"third"` | 0.0442 | first cousins, and the above |

Pairs with a missing kinship estimate are treated as unrelated, with a
warning.

The sample check runs here, rather than in the selectors, because it
needs the kinship of every pair of samples, including the unrelated
ones, which the network no longer holds. A network built in another way
can be checked by calling
[`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md)
on the kinship estimates it was built from.

## References

Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010).
Robust relationship inference in genome-wide association studies.
*Bioinformatics*, 26(22), 2867–2873.
[doi:10.1093/bioinformatics/btq559](https://doi.org/10.1093/bioinformatics/btq559)

## See also

[`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md),
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md),
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md),
[`rkn_select_global`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_global.md),
[`rkn_plot_pre`](https://paolodalena.github.io/RedKinNet/reference/rkn_plot_pre.md).

## Examples

``` r
# Kinship estimates for six samples: a parent and child (A, B), first cousins (C, D)
# and two unrelated samples (E, F).
K <- matrix(-0.01, 6, 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
K["A", "B"] <- K["B", "A"] <- 0.25
K["C", "D"] <- K["D", "C"] <- 0.06
diag(K) <- 0.5

# Relatives of the third degree or closer: two families of two, two singletons
g <- rkn_graph_from_kinship(K, threshold = "third")
igraph::components(g)$csize
#> [1] 2 2 1 1

# First-degree relatives only: the cousins count as unrelated
igraph::components(rkn_graph_from_kinship(K, threshold = "first"))$csize
#> [1] 2 1 1 1 1

# The same data as a pair table that lists only the related pairs, as KING does;
# samples adds the unrelated E and F
pairs <- data.frame(ID1 = c("A", "C"), ID2 = c("B", "D"), Kinship = c(0.25, 0.06))
g <- rkn_graph_from_kinship(pairs, threshold = "third", samples = LETTERS[1:6])

# Covariates for a criteria-based selection
covariates <- data.frame(id = LETTERS[1:6], sex = c("M", "F", "F", "M", "M", "F"))
g <- rkn_graph_from_kinship(K, threshold = "third", covariates = covariates)
rkn_select_criteria(g, rkn_criterion_level("sex", "F"), seed = 1)$kept
#> [1] "B" "C" "E" "F"
```
