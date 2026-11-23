# Plot scaled compositions

Sums the draws from
[`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md)
over all strata within levels of `grps` and plots the resulting
compositions (as proportions or as catch-at-bin) with uncertainty
intervals. Optional reference compositions (e.g. raw samples or a known
truth) can be overlaid.

## Usage

``` r
scaled_comp_plot(
  scaled_comp = NULL,
  grps = NULL,
  scaled = T,
  scales = "free_y",
  comp_are = NULL,
  comp_are2 = NULL,
  cvar = "tcatch",
  xlab = "Length bin (cm)",
  comp_label = "Reference",
  comp2_label = "Reference 2"
)
```

## Arguments

- scaled_comp:

  output of
  [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md)

- grps:

  grouping variables to aggregate to (e.g. `"Year"`)

- scaled:

  If `TRUE`, plot proportions; if `FALSE`, plot catch-at-bin.

- scales:

  Passed to `facet_wrap()`/`facet_grid()`. Panels are a `facet_wrap()`
  by `grps` when it has one element, and a `facet_grid()` with the first
  element in rows and the remaining elements in columns otherwise.

- comp_are, comp_are2:

  Optional data frames with `grps`, `bin` and `cvar`, overlaid as blue
  and orange points respectively. Counts are converted to proportions
  within `grps` when `scaled = TRUE`.

- cvar:

  Column in `comp_are`/`comp_are2` holding the counts or catch per bin.

- xlab:

  x-axis label (the bin variable)

- comp_label, comp2_label:

  Legend labels for the points from `comp_are` (skyblue) and `comp_are2`
  (orange).

## Value

A ggplot with the mean and 95% quantile interval of the composition (or
catch-at-bin) per bin.

## Examples

``` r
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                 iter = 300, warmup = 100, thin = 2)
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)
sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))

scaled_comp_plot(sc, grps = "Year")


# overlay the raw sample composition
raw <- dplyr::rename(ex_comp, bin = length)
scaled_comp_plot(sc, grps = "Year", comp_are = raw, cvar = "count")
```
