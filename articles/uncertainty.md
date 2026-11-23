# Selection uncertainty

A cohort usually has an enormous number of maximum sets, all of the same
size. Choosing one of them is arbitrary, and the individuals analysed,
and therefore the results of the analysis, depend on that choice. This
vignette shows how RedKinNet quantifies this uncertainty on the example
cohort, with the network built as in *Get started with RedKinNet*.

``` r

library(RedKinNet)
kin <- rkn_example$kinship
covariates <- rkn_example$samples
chk <- rkn_check_samples(kin, threshold = "third", samples = covariates$id)
g <- rkn_graph_from_kinship(kin, threshold = "third", samples = covariates$id,
                            exclude = chk$contamination$id, covariates = covariates,
                            check = FALSE)
```

## How many maximum sets

A uniform selection counts the maximum sets of every family exactly. The
number for the whole network is the product of the numbers of the
families.

``` r

sel <- rkn_select_random(g, uniform = TRUE, seed = 1)
sel
#> <rkn_selection> Uniform random selection, seed 1
#> Call: rkn_select_random(g = g, uniform = TRUE, seed = 1)
#> Individuals: 751 | selected: 452 (60.2%)
#> Families: 376 | size min/mean/max: 1 / 2.00 / 31
#> Families with up to 5 members: 357 | with 6 or more: 19
#> Families sampled uniformly among their maximum sets: 376 of 376
#> Maximum sets: 4.4e+54 in total; families with more than one: 127 of 376
#> Two uniform draws differ on average in 84.6 of the 452 selected individuals (18.7%)
fs <- sel$family_summary
head(fs[order(-fs$n_max_sets), c("family_id", "family_size", "n_kept", "n_max_sets")])
#>    family_id family_size n_kept n_max_sets
#> 21        21          16      5         40
#> 10        10          10      3         18
#> 18        18          31      6         16
#> 38        38          10      4         13
#> 41        41           7      2          9
#> 57        57           6      2          7
```

## Inclusion probabilities

The inclusion probability of an individual is the proportion of the
maximum sets that contain it, which is its probability of being selected
by a uniform draw. Individuals with probability 1 are selected in every
maximum set, those with probability 0 in none; only the others are
actually chosen.

``` r

p <- sel$node_diagnostics$inclusion_probability
table(cut(p, c(-0.01, 0, 0.25, 0.5, 0.75, 0.99, 1),
          labels = c("0", "(0, 0.25]", "(0.25, 0.5]", "(0.5, 0.75]", "(0.75, 1)", "1")))
#> 
#>           0   (0, 0.25] (0.25, 0.5] (0.5, 0.75]   (0.75, 1)           1 
#>          79          78         274          21           2         297
```

Two independent uniform draws differ, on average, in the sum of
`p * (1 - p)` selected individuals, the figure reported by
[`print()`](https://rdrr.io/r/base/print.html).

``` r

sum(p * (1 - p))
#> [1] 84.55876
```

## How much a result can change

Repeating the uniform draw with different seeds shows how much a
statistic computed on the selected individuals varies across equally
valid selections. Here the statistics are the proportion of women and
the mean age; in a real study they would be the results of the
downstream analysis. Each seed gives an independent draw.

``` r

draws <- data.frame(seed = 1:50, women = NA_real_, mean_age = NA_real_)
for (s in draws$seed) {
  kept <- rkn_select_random(g, uniform = TRUE, seed = s, diagnostics = FALSE)$kept
  i <- match(kept, covariates$id)
  draws$women[s] <- mean(covariates$sex[i] == "F", na.rm = TRUE)
  draws$mean_age[s] <- mean(covariates$age[i], na.rm = TRUE)
}
summary(draws[, c("women", "mean_age")])
#>      women           mean_age    
#>  Min.   :0.4855   Min.   :49.78  
#>  1st Qu.:0.5000   1st Qu.:50.26  
#>  Median :0.5056   Median :50.54  
#>  Mean   :0.5049   Mean   :50.57  
#>  3rd Qu.:0.5101   3rd Qu.:50.87  
#>  Max.   :0.5237   Max.   :51.69
```

The spread of these values is part of the uncertainty of any result
obtained on one selection. Reporting the seed, the number of maximum
sets and the average difference between draws makes the selection
reproducible and its arbitrariness visible; repeating the analysis over
several draws shows whether a conclusion depends on it.

## Why uniform

With `uniform = FALSE`,
[`rkn_select_random()`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
resolves larger families by structural rules, which never choose some of
their maximum sets. On a chain of six relatives, for instance, there are
four maximum sets, and the rules choose only between the two of
alternate relatives. A uniform draw gives every maximum set the same
chance, so that the inclusion probabilities above describe the selection
exactly.

``` r

chain <- igraph::make_graph(~ A - B - C - D - E - F)
rules <- character(200)
uniform <- character(200)
for (s in 1:200) {
  rules[s] <- paste(sort(rkn_select_random(chain, uniform = FALSE, seed = s, diagnostics = FALSE)$kept),
                    collapse = "")
  uniform[s] <- paste(sort(rkn_select_random(chain, uniform = TRUE, seed = s, diagnostics = FALSE)$kept),
                      collapse = "")
}
table(rules)
#> rules
#> ACE BDF 
#>  95 105
table(uniform)
#> uniform
#> ACE ACF ADF BDF 
#>  55  49  46  50
```
