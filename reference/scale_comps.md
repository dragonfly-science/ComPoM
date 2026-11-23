# Scale predicted compositions to catch

Predicts the composition for each stratum in `scale_df` from a fitted
model and multiplies the predicted proportions by the stratum weight
(`predvar`, e.g. catch or predicted CPUE). Strata are the unique
combinations of `grps` (plus the coordinates for spatial models). The
result is a set of draws of numbers-at-bin per stratum, which can be
summed over any subset of `grps` and normalised to obtain weighted
compositions (e.g. with
[`scaled_comp_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_comp_plot.md)).

## Usage

``` r
scale_comps(
  scale_df,
  predvar = "catch",
  fit = NULL,
  grps = NULL,
  iters = NULL,
  pgrid = NULL,
  remove_nuisance = FALSE
)
```

## Arguments

- scale_df:

  data frame used for scaling, with one or more rows per stratum,
  containing `grps`, `predvar`, any `cvars`, and (spatial fits) the
  coordinate and time columns.

- predvar:

  Column in `scale_df` used for scaling (e.g. catch)

- fit:

  A `ComPoM` object from
  [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)

- grps:

  Grouping variables defining the scaling strata. Must include the
  `ffx_form` variables and, for catch weighting, the `re_form`
  variables.

- iters:

  Number of posterior draws (brms only; default uses all stored draws).
  For sdmTMB fits the number of draws is set by
  [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md).

- pgrid:

  Not used; reserved for future use.

- remove_nuisance:

  Logical. `FALSE` (default) selects catch weighting, with nuisance
  random effects evaluated at each stratum's level. `TRUE` selects
  CPUE-based scaling, with nuisance random effects removed.

## Value

A data frame of draws (`.draw`) with predicted numbers-at-bin
(`tot_by_bin`) per stratum and bin. The weight per stratum (sum of
`predvar`) is in `n`.

## Details

The two settings of `remove_nuisance` correspond to two uses:

- Catch weighting (`remove_nuisance = FALSE`, default):

  The aim is the composition of the total catch of the full fleet. All
  model terms are evaluated at the values of each stratum in `scale_df`,
  including the nuisance random effects (`re_form`). Every nuisance
  variable must therefore be in `grps`. `predvar` is the catch.

- CPUE-based scaling (`remove_nuisance = TRUE`):

  The aim is a composition representative of the population,
  standardised for variables that affect the local availability of sizes
  to the sampling (e.g. fleet, season). Nuisance random effects are
  removed (population level); fixed effects, smooths and spatial fields
  are retained. `predvar` is typically predicted CPUE or abundance per
  stratum.

In both cases:

- Fixed effects (`ffx_form`) are evaluated at the levels in `scale_df`;
  these variables must therefore be in `grps`.

- Continuous covariates (`cvars`) are evaluated at their mean within
  each stratum of `scale_df`. Setting them to a reference value in
  `scale_df` standardises the compositions to that value.

- Spatial and spatio-temporal fields (spatial sdmTMB fits) are evaluated
  at the coordinates in `scale_df`.

Nuisance levels that were not in the fitted data (e.g. unsampled
vessels) are handled differently by the backends. brms uses, in each
draw, the posterior draw of a randomly chosen fitted level
(`sample_new_levels = "uncertainty"`). sdmTMB predicts such rows at the
population level of all nuisance random effects, and a warning is
issued.

## See also

[`scaled_comp_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_comp_plot.md),
[`scaled_ridge_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_ridge_plot.md)

## Examples

``` r
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                 iter = 300, warmup = 100, thin = 2)
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)

# catch weighting: composition of the total catch of the fleet
sc <- scale_comps(ex_catch, predvar = "catch", fit = fit,
                  grps = c("Year", "Season", "Vessel"))
head(sc)
#> # A tibble: 6 × 9
#>   Year  Season Vessel bin       n  .row iter  tot_by_bin .draw
#>   <fct> <fct>  <fct>  <fct> <dbl> <int> <chr>      <dbl> <dbl>
#> 1 2015  summer V1     10     1904     1 V1          436.     1
#> 2 2015  summer V1     10     1904     1 V2          528.     2
#> 3 2015  summer V1     10     1904     1 V3          357.     3
#> 4 2015  summer V1     10     1904     1 V4          498.     4
#> 5 2015  summer V1     10     1904     1 V5          412.     5
#> 6 2015  summer V1     10     1904     1 V6          436.     6
scaled_comp_plot(sc, grps = "Year")


# CPUE-based scaling: nuisance effects removed; one weight per year
cpue <- data.frame(Year = factor(2015:2020), cpue = 1)
sc_pop <- scale_comps(cpue, predvar = "cpue", fit = fit, grps = "Year",
                      remove_nuisance = TRUE)
```
