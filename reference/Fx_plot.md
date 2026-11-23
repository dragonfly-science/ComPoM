# Plot bin-specific effects as multipliers

Predicts the composition for each level of `grp`, holding all other
covariates at their mean (numeric) or mode (factor), and with the offset
fixed. Within each level and draw, predictions are converted to
proportions. Each proportion is then divided by the geometric mean of
that bin's proportions across levels of `grp` and draws. A multiplier
above 1 indicates that the bin is relatively more frequent in that level
than on average.

## Usage

``` r
Fx_plot(
  mod,
  grp = "Year",
  form = "~(1|bin:Year)",
  grid = NULL,
  cvar = NULL,
  xlab = "Length bin (cm)"
)
```

## Arguments

- mod:

  An object from fit_model

- grp:

  Character vector of the variable(s) to plot effects for (e.g.
  `"Year"`); must include `cvar` if one is used. Facets are by the first
  variable.

- form:

  Unused; retained for backwards compatibility.

- grid:

  Optional character vector of two variables in `grp` for a
  `facet_grid()` layout (rows, columns) instead of `facet_wrap()`.

- cvar:

  Name of the continuous covariate, if the model has one. It is held at
  its mean.

- xlab:

  x-axis label (the bin variable)

## Value

A ggplot of median multipliers and intervals per bin, faceted by the
first variable in `grp` (or by `grid`).

## Examples

``` r
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                 iter = 300, warmup = 100, thin = 2)
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)
Fx_plot(fit, grp = "Year")
```
