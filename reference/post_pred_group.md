# Posterior predictive check by group

Compares observed counts per bin with posterior predictive (brms) or
simulated (sdmTMB) counts, summed within levels of `grp`. Predicted
totals are rescaled to the observed total within each group and draw, so
the check is on the composition rather than on the total. Requires a fit
from `fit_model(..., add_preds = TRUE)`.

## Usage

``` r
post_pred_group(mod, grp = NULL, xlab = "Length (cm)")
```

## Arguments

- mod:

  An object from fit_model

- grp:

  Name(s) of the grouping variable(s) in the model data to check over
  (e.g. `"Year"`). `NULL` pools all strata.

- xlab:

  x-axis label (the bin variable)

## Value

A ggplot with predicted medians and intervals (black) and observed
counts (blue), faceted by `grp`.

## See also

[`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)

## Examples

``` r
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                 iter = 300, warmup = 100, thin = 2, add_preds = TRUE)
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)
post_pred_group(fit, grp = "Year", xlab = "Length (cm)")
```
