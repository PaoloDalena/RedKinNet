# Structural rules for families of up to 5 members

A curated list of structural rules used by
[`rkn_select_random`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_random.md)
and
[`rkn_select_criteria`](https://paolodalena.github.io/RedKinNet/reference/rkn_select_criteria.md)
on families of up to 5 members: each rule recognises a family structure,
keeps or discards some members outright, and groups the others into the
alternatives left to choose from. Every possible structure of a family
of up to 5 members matches exactly one rule. The rules are not used in
uniform selections (`rkn_select_random(uniform = TRUE)`) nor in
selections over all maximum sets
(`rkn_select_criteria(optimal = TRUE)`).

## Usage

``` r
rules_db
```

## Format

A list of rule objects (lists).

## Details

Each element of the list is a rule with fields such as:

- `name`

- `size_allowed`

- `match_fn`

- `kept_rules`

- `discarded_rules`

- `conflict_rules`

- optional `use_largest_ivs` (logical): if `TRUE`, the remaining
  ambiguity is resolved by choosing among all the maximum sets of the
  eligible members

- optional `kept_fn`, `discarded_fn`

## Examples

``` r
length(rules_db)
#> [1] 30
rules_db[[2]]$name
#> [1] "pair_edge"
```
