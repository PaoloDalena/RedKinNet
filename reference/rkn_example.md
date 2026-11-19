# Example cohort with simulated kinship estimates

A synthetic cohort of 752 samples for the examples and vignettes of the
package: a table of pairwise kinship estimates, in the format written by
KING or by PLINK 2 when the output is restricted to pairs with a kinship
of at least 0.02, and a table of covariates. The cohort was simulated
from pedigrees: independent families of up to four generations, of which
few members were sampled, and a small isolated village whose sampled
members are related through several lines of descent. One sample was
genotyped twice and one sample is contaminated;
[`rkn_check_samples`](https://paolodalena.github.io/RedKinNet/reference/rkn_check_samples.md)
finds both. No real individual is involved.

## Usage

``` r
rkn_example
```

## Format

A list of two data frames:

- `kinship`:

  the pairs of samples with an estimated kinship of at least 0.02, with
  columns `id1`, `id2` and `kinship`. Samples without relatives appear
  in no row, so the identifiers of all the samples are passed
  separately, as `rkn_example$samples$id`.

- `samples`:

  one row per sample, with `id`; `sex`, `"F"` or `"M"`; `age`, in years;
  and `call_rate`, the proportion of genotypes called. Sex and age have
  a few missing values.

## Source

Simulated with the script `data-raw/rkn_example.R` of the package
sources: pedigree kinship with variation of the realised relatedness and
estimation noise.

## Examples

``` r
str(rkn_example$samples)
#> 'data.frame':    752 obs. of  4 variables:
#>  $ id       : chr  "S001" "S002" "S003" "S004" ...
#>  $ sex      : chr  NA "M" "M" "F" ...
#>  $ age      : int  61 21 48 86 60 58 87 58 57 57 ...
#>  $ call_rate: num  0.998 0.981 0.999 0.998 0.999 ...
head(rkn_example$kinship)
#>    id1  id2 kinship
#> 1 S001 S047  0.0369
#> 2 S001 S055  0.1139
#> 3 S001 S578  0.1364
#> 4 S001 S620  0.2500
#> 5 S001 S672  0.0654
#> 6 S002 S323  0.0557

# Check the samples, then build the network without the contaminated sample
chk <- rkn_check_samples(rkn_example$kinship, threshold = "third",
                         samples = rkn_example$samples$id)
chk
#> <rkn_sample_check> 752 samples, relatedness threshold 0.0442
#> 
#> Possible duplicates or monozygotic twins (kinship >= 0.354): 1 pair(s)
#>   id1  id2 kinship
#>  S685 S707   0.496
#> 
#> Possible contamination: 1 sample(s)
#> (related to unusually many samples)
#>    id n_relatives mean_kinship z_relatives z_kinship largest_family_without
#>  S672         631           NA         8.9        NA                     31
#> Largest family of the network: 717 samples;
#> without all the flagged samples together: 31 samples.
#> 
#> Ancestry not assessed: the mean kinship needs every pair of samples.
g <- rkn_graph_from_kinship(rkn_example$kinship, threshold = "third",
                            samples = rkn_example$samples$id,
                            exclude = chk$contamination$id,
                            covariates = rkn_example$samples, check = FALSE)
rkn_select_random(g, uniform = TRUE, seed = 1)
#> <rkn_selection> Uniform random selection, seed 1
#> Call: rkn_select_random(g = g, uniform = TRUE, seed = 1)
#> Individuals: 751 | selected: 452 (60.2%)
#> Families: 376 | size min/mean/max: 1 / 2.00 / 31
#> Families with up to 5 members: 357 | with 6 or more: 19
#> Families sampled uniformly among their maximum sets: 376 of 376
#> Maximum sets: 4.4e+54 in total; families with more than one: 127 of 376
#> Two uniform draws differ on average in 84.6 of the 452 selected individuals (18.7%)
```
