# Prepare compositional count data

Takes a dataset of counts (number of sampled individuals) per bin
(length/age/etc) across a range of strata (`vars_for_grouping`) and
expands it for use in the Poisson-multinomial model. Counts are summed
within stratum and bin, bins without observations within a stratum are
filled with zero counts, and the stratum total `n` (used as offset in
[`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md))
is added.

## Usage

``` r
data_prep(
  comp = NULL,
  vars_for_grouping = NULL,
  bin_lab = "bin",
  count_var = "count",
  min_n = 1L,
  bin_width = NULL
)
```

## Arguments

- comp:

  A data frame with one row per observation (or per stratum and bin),
  with a column of counts, a column of bins and all `vars_for_grouping`.

- vars_for_grouping:

  Character vector of the variables that define a stratum (e.g.
  `c("Year", "Season", "Vessel")`).

- bin_lab:

  Name of the column holding the bins (length, age, etc.). Default
  `"bin"`. The column is renamed to `bin` in the output.

- count_var:

  Name of the column holding the counts. Default `"count"`.

- min_n:

  Integer. Minimum total count per grouping stratum. Strata with fewer
  than `min_n` total fish are dropped before fitting. Useful for
  removing sparse early years that would create all-zero design matrix
  columns and a singular Hessian. Default `1L` (only empty strata
  dropped).

- bin_width:

  Numeric. If supplied, bins are aggregated into wider bins of this
  width before summarising (e.g. `bin_width = 20` merges 10cm bins into
  20cm bins). Applied before `min_n` filtering. Default `NULL` (no
  aggregation).

## Value

A grouped data frame with one row per stratum and bin, containing
`vars_for_grouping`, `bin`, the count per bin (`tot_by_bin`) and the
stratum total (`n`).

## Details

Strata are defined by the unique combinations of `vars_for_grouping`.
Any variable later used in
[`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
(fixed effects, nuisance random effects, continuous covariates,
coordinates, time) must be included here, otherwise it is aggregated
away.

## See also

[`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)

## Examples

``` r
dat <- data_prep(
  ex_comp,
  vars_for_grouping = c("Year", "Season", "Vessel", "trip"),
  bin_lab = "length",
  count_var = "count"
)
head(dat)
#> # A tibble: 6 × 7
#> # Groups:   Year, Season, Vessel, trip [1]
#>   Year  Season Vessel  trip   bin tot_by_bin     n
#>   <fct> <fct>  <fct>  <int> <dbl>      <int> <int>
#> 1 2015  summer V1         1    10          5    60
#> 2 2015  summer V1         1    20         15    60
#> 3 2015  summer V1         1    30         21    60
#> 4 2015  summer V1         1    40         15    60
#> 5 2015  summer V1         1    50          2    60
#> 6 2015  summer V1         1    60          2    60

# merge 10 cm bins into 20 cm bins
head(data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count",
               bin_width = 20))
#> # A tibble: 6 × 7
#> # Groups:   Year, Season, Vessel, trip [2]
#>   Year  Season Vessel  trip   bin tot_by_bin     n
#>   <fct> <fct>  <fct>  <int> <dbl>      <int> <int>
#> 1 2015  summer V1         1     0          5    60
#> 2 2015  summer V1         1    20         36    60
#> 3 2015  summer V1         1    40         17    60
#> 4 2015  summer V1         1    60          2    60
#> 5 2015  summer V1         2     0         21    60
#> 6 2015  summer V1         2    20         33    60
```
