# Influence of model terms on the composition

For each term, computes (1) the effect of each level on every bin, (2)
the distribution of the data over the levels by year, and (3) the
influence of the term on the annual composition. The influence is the
log ratio of the annual proportion in each bin predicted with the term
at the observed levels to the proportion predicted with the term at its
reference. The same quantity is also given for a statistic of the
composition (by default the mean bin). The design follows the
coefficient-distribution-influence plots used for CPUE standardisation.

## Usage

``` r
influ_comps(
  fit,
  scale_df = NULL,
  predvar = "catch",
  grps,
  by,
  terms = NULL,
  grid = NULL,
  selectivity = NULL,
  reference = NULL,
  weight = NULL,
  density = NULL,
  stat = "mean",
  components = c("catch", "index")
)
```

## Arguments

- fit:

  A `ComPoM` object from
  [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
  with `backend = "TMB"`.

- scale_df, predvar:

  Data frame and weight column for the catch component; see
  [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md).

- grps:

  Grouping variables for the scaling; see
  [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md).

- by:

  Variable in `grps` (or `grid`) over which the statistic is plotted,
  usually the year.

- terms:

  Terms to evaluate. Default: the `ffx_form` and `re_form` variables
  except `by`.

- grid, selectivity, reference, weight, density:

  Settings of the index component; see
  [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md).

- stat:

  `"mean"` (default), a probability between 0 and 1 for a quantile of
  the bin (e.g. `0.9`), or a function of `bin` and `p`. See
  [`comp_quantile()`](https://dragonfly-science.github.io/ComPoM/reference/comp_quantile.md).

- components:

  Components to compute: `"catch"` and/or `"index"`.

## Value

A list of class `comp_influ` with data frames `effects` (log multiplier
by term, level and bin), `distribution` (share of the weight by term,
level, component and year), `influence` (log ratio by term, component,
year and bin) and `influence_stat` (difference in the statistic by term,
component and year).

## Details

The influence is computed for the catch component (`mode = "catch"` in
[`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md))
and the index component (`mode = "index"`). Terms in `selectivity` have
no influence on the index by construction and are given zero influence
for that component. Random effects are removed with the TMB backend only
together, so that with more than one random-effect variable they form
one term (labelled with all variables). Effects of covariate smooths and
spatial fields are not included; see
[`cvar_smooth_plot()`](https://dragonfly-science.github.io/ComPoM/reference/cvar_smooth_plot.md)
and
[`spatial_field_plot()`](https://dragonfly-science.github.io/ComPoM/reference/spatial_field_plot.md).

## See also

[`plot_cdi_comps()`](https://dragonfly-science.github.io/ComPoM/reference/plot_cdi_comps.md),
[`step_comps()`](https://dragonfly-science.github.io/ComPoM/reference/step_comps.md)

## Examples

``` r
# \donttest{
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
fit <- fit_model("Year + Season", "Vessel", data = dat, backend = "TMB",
                 iter = 300, warmup = 100, thin = 2)
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + factor(Season):bin + (1 | bin:Vessel)
inf <- influ_comps(fit, scale_df = ex_catch, predvar = "catch",
                   grps = c("Year", "Season", "Vessel"), by = "Year",
                   terms = "Season")
plot_cdi_comps(inf, "Season")

# }
```
