# Check samples for duplicates, possible contamination and different ancestry

Summarises, for each sample, its relatives and its mean kinship with the
rest of the cohort, and flags three kinds of anomalies before a
selection:

- **possible duplicates or monozygotic twins**: pairs with a kinship of
  at least 0.354. They need no action for the selection, which keeps at
  most one of them, but they are worth a quality check;

- **possible contamination**: samples related to unusually many others,
  with an unusually high kinship even with the samples they are not
  related to. A sample containing DNA from more than one person tends to
  appear slightly related to everyone, and can merge many families into
  one; excluding it (argument `exclude` of
  [`rkn_graph_from_kinship`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md))
  restores the real families;

- **possible different ancestry**: samples with an unusually low kinship
  with the samples they are not related to. The ancestry of an
  individual is the population their forebears come from; populations
  separated for many generations differ in the frequencies of genetic
  variants, and KING estimates are negative for pairs of individuals
  from different populations (Manichaikul et al., 2010). This flag is
  informative only: these samples are related to no one, so they change
  neither the network nor the selection, and they need not be excluded
  for it. It is reported because ancestry matters for some downstream
  analyses, such as association studies, where it is usually accounted
  for with the principal components of the genotypes (Price et al.,
  2006). A low mean kinship is only a hint: assessing ancestry requires
  the genotypes.

Nothing is removed: the decision is left to the user.
[`rkn_graph_from_kinship`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md)
runs this check by default (argument `check`) and summarises the
duplicates and the contamination it finds in a warning.

## Usage

``` r
rkn_check_samples(
  x,
  threshold,
  samples = NULL,
  z_kinship = 5,
  z_relatives = 3,
  z_relatives_only = 5,
  z_ancestry = -5
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

- z_kinship, z_relatives:

  A sample is flagged as possible contamination when the robust z-score
  of its mean kinship exceeds `z_kinship` and that of its number of
  relatives exceeds `z_relatives`.

- z_relatives_only:

  When the mean kinship cannot be computed, a sample is flagged as
  possible contamination when the robust z-score of its number of
  relatives exceeds `z_relatives_only`.

- z_ancestry:

  A sample is flagged as possible different ancestry when the robust
  z-score of its mean kinship is below `z_ancestry`.

## Value

An object of class `rkn_sample_check`, with a `print` method: a list
with elements

- `samples`:

  a data frame with one row per sample: `id`; `n_relatives`, the number
  of samples with a kinship of at least `threshold`; `n_duplicate`,
  `n_first`, `n_second` and `n_third`, the relatives by degree;
  `mean_kinship`, the mean kinship with the samples it is not related to
  (see Details); the robust z-scores `z_relatives` and `z_kinship`;
  `flag`, `"possible contamination"`, `"possible different ancestry"` or
  an empty string;

- `duplicates`:

  a data frame of the pairs with a kinship of at least 0.354;

- `contamination`:

  a data frame of the samples flagged as possible contamination, with
  `largest_family_without`, the size of the largest family of the
  network after removing the sample;

- `threshold`, `largest_family`:

  the threshold used and the size of the largest family of the network;

- `largest_family_without_flagged`:

  the size of the largest family after removing all the samples flagged
  as possible contamination together. When several contaminated samples
  are related to the same samples, removing one of them leaves the
  others holding the families together, so only the joint removal shows
  their effect. It equals `largest_family` when no sample is flagged;

- `mean_kinship_available`:

  whether the mean kinship could be computed.

## Details

The robust z-score of a value is its distance from the median of the
cohort in units of the median absolute deviation (MAD), so that the
anomalies being sought do not inflate the scale. When the MAD is zero,
as for the number of relatives in a cohort where most samples have none,
the mean absolute deviation from the median is used instead, rescaled to
be comparable with the MAD. The number of relatives enters on the
logarithmic scale, `log(1 + relatives)`.

The mean kinship of a sample is computed over the samples it is not
related to, that is over the pairs below the third-degree cutoff of KING
(0.0442): true relatives raise the kinship of a few pairs only, whereas
contamination raises it with everyone, so leaving the relatives out
keeps a sample with many true relatives, or a member of a small cohort
with one close relative, from looking anomalous. The mean therefore
requires a kinship matrix or a pair table that lists every pair. With a
table of related pairs only, it is not available: contamination is then
judged on the number of relatives alone (see `z_relatives_only`) and
ancestry is not assessed.

Relatives are counted by degree with the KING cutoffs listed in
[`rkn_graph_from_kinship`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md).

## References

Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010).
Robust relationship inference in genome-wide association studies.
*Bioinformatics*, 26(22), 2867–2873.
[doi:10.1093/bioinformatics/btq559](https://doi.org/10.1093/bioinformatics/btq559)

Price AL, Patterson NJ, Plenge RM, Weinblatt ME, Shadick NA, Reich D
(2006). Principal components analysis corrects for stratification in
genome-wide association studies. *Nature Genetics*, 38(8), 904–909.
[doi:10.1038/ng1847](https://doi.org/10.1038/ng1847)

## See also

[`rkn_graph_from_kinship`](https://paolodalena.github.io/RedKinNet/reference/rkn_graph_from_kinship.md)

## Examples

``` r
# A synthetic cohort of 200 samples: 30 parent-child pairs, one sample genotyped
# twice (S199 and S200), one contaminated sample (S001), whose kinship is raised with
# everyone, and 15 samples of a different ancestry (S181 to S195).
set.seed(1)
ids <- sprintf("S%03d", 1:200)
K <- matrix(rnorm(200^2, -0.015, 0.004), 200, 200, dimnames = list(ids, ids))
K[upper.tri(K)] <- t(K)[upper.tri(K)]
ancestry <- 181:195
K[ancestry, -ancestry] <- rnorm(15 * 185, -0.15, 0.03)
K[-ancestry, ancestry] <- t(K[ancestry, -ancestry])
for (i in seq(2, 60, by = 2)) K[i, i + 1] <- K[i + 1, i] <- 0.25
K["S199", "S200"] <- K["S200", "S199"] <- 0.49
K["S001", -1] <- K[-1, "S001"] <- rnorm(199, 0.045, 0.01)
diag(K) <- 0.5

chk <- rkn_check_samples(K, threshold = "third")
chk
#> <rkn_sample_check> 200 samples, relatedness threshold 0.0442
#> 
#> Possible duplicates or monozygotic twins (kinship >= 0.354): 1 pair(s)
#>   id1  id2 kinship
#>  S199 S200    0.49
#> 
#> Possible contamination: 1 sample(s)
#> (related to unusually many samples,
#> with an unusually high kinship even with unrelated ones)
#>    id n_relatives mean_kinship z_relatives z_kinship largest_family_without
#>  S001         105       0.0364        10.2      81.5                      2
#> Largest family of the network: 121 samples;
#> without all the flagged samples together: 2 samples.
#> 
#> Possible different ancestry (mean kinship far below the cohort): 15 sample(s)
#> These do not affect the selection; see $samples for the values.
head(chk$samples)
#>     id n_relatives n_duplicate n_first n_second n_third mean_kinship
#> 1 S001         105           0       0        0     105   0.03636107
#> 2 S002           1           0       1        0       0  -0.02538556
#> 3 S003           2           0       1        0       1  -0.02546515
#> 4 S004           2           0       1        0       1  -0.02614858
#> 5 S005           1           0       1        0       0  -0.02548009
#> 6 S006           2           0       1        0       1  -0.02459609
#>   z_relatives  z_kinship                   flag
#> 1   10.152469 81.5466443 possible contamination
#> 2    0.000000 -0.3442226                       
#> 3    1.036818 -0.4497839                       
#> 4    1.036818 -1.3561702                       
#> 5    0.000000 -0.4695907                       
#> 6    1.036818  0.7028056                       

# Build the network without the possibly contaminated sample; the warning still
# reports the duplicates
g <- rkn_graph_from_kinship(K, threshold = "third", exclude = chk$contamination$id)
#> Warning: Sample check: 1 pair of possible duplicates or monozygotic twins. See rkn_check_samples() for details and the argument `exclude` to leave samples out; check = FALSE skips this check.
```
