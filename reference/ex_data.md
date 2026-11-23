# Simulated example data

Simulated length-composition sampling of a fishery, used in the examples
and the README. Fish get larger over the years, vessels select different
sizes, winter fish are larger, fish size varies along a spatial
gradient, and sea-surface temperature shifts the composition. `ex_comp`
holds the sampled length frequencies and `ex_catch` the catch of the
strata, for use with
[`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md).

## Format

`ex_comp` is a data frame with one row per stratum and length bin (144
strata x 6 bins):

- Year:

  Factor, 2015 to 2020.

- Season:

  Factor, `"summer"` or `"winter"`.

- Vessel:

  Factor, `"V1"` to `"V4"`.

- trip:

  Trip number within year, season and vessel.

- x, y:

  Sampling location (arbitrary projected coordinates, 0 to 100).

- sst:

  Sea-surface temperature at the sample.

- length:

  Length bin (10 to 60).

- count:

  Number of fish measured in the bin (60 fish per stratum).

`ex_catch` is a data frame with one row per stratum (144 rows) and the
stratum variables of `ex_comp` (`Year`, `Season`, `Vessel`, `trip`, `x`,
`y`, `sst`) plus `catch`, the total catch of the stratum.

## Source

Simulated, see `data-raw/ex_data.R`.

## Examples

``` r
head(ex_comp)
#>   Year Season Vessel trip        x        y   sst length count
#> 1 2015 summer     V1    1 86.37275 36.07672 14.02     10     5
#> 2 2015 summer     V1    1 86.37275 36.07672 14.02     20    15
#> 3 2015 summer     V1    1 86.37275 36.07672 14.02     30    21
#> 4 2015 summer     V1    1 86.37275 36.07672 14.02     40    15
#> 5 2015 summer     V1    1 86.37275 36.07672 14.02     50     2
#> 6 2015 summer     V1    1 86.37275 36.07672 14.02     60     2
head(ex_catch)
#>   Year Season Vessel trip         x         y   sst catch
#> 1 2015 summer     V1    1 86.372746 36.076716 14.02   250
#> 2 2016 summer     V1    1  6.476116 14.744300 15.94   508
#> 3 2017 summer     V1    1 86.953293 41.618799 19.06   534
#> 4 2018 summer     V1    1 13.818807 47.281130 19.46   581
#> 5 2019 summer     V1    1 75.313348 18.731764 19.26   692
#> 6 2020 summer     V1    1 22.750807  4.632714 17.18   370
```
