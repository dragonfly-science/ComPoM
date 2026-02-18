Install the development version from source:

```r
# install.packages("devtools")
devtools::install_github("dragonfly-science/ComPoM")
```

Or from a local clone:

```r
devtools::install()lation

You can install the development version from your local source using:

```r
# In the package root directory
devtools::install()
```

Or, if available on GitHub:

```r
# Replace 'yourusername/ComPoM' with the actual repo
devtools::install_github("yourusername/ComPoM")
```

## Capabilities

| Function | Description |
|---|---|
| `data_prep()` | Prepares compositional count data for modelling |
| `parse_form()` | Constructs the Poisson-multinomial model formula |
| `fit_model()` | Fits the model using `brms` or `sdmTMB` backends |
| `Fx_plot()` | Plots compositional effects by group |
| `scale_comps()` | Scales compositions by a catch or effort variable |
| `scaled_comp_plot()` | Plots scaled compositions with uncertainty |
| `scaled_ridge_plot()` | Ridgeline plot of scaled compositions |
| `post_pred_group()` | Posterior predictive check by group |

## Basic Usage

```r
library(ComPoM)

# 1. Prepare compositional count data
#    comp: data frame with columns for bin, count, and grouping variables
dat <- data_prep(
  comp             = my_data,
  vars_for_grouping = c("year", "area"),
  bin_lab          = "length_bin",
  count_var        = "count"
)

# 2. Fit the model (brms backend by default)
mod <- fit_model(
  form    = "year + area",
  data    = dat,
  backend = "brms",
  chains  = 4,
  iter    = 2000
)

# 3. Plot compositional effects
Fx_plot(mod, grp = "year")

# 4. Scale compositions by catch and plot
scaled <- scale_comps(
  scale_df = catch_data,
  predvar  = "catch",
  fit      = mod,
  grps     = c("year", "area")
)

scaled_comp_plot(scaled, grps = c("year", "area"))
```

### Spatial model (sdmTMB backend)

```r
library(sdmTMB)

mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 20)

mod_spatial <- fit_model(
  form    = "year + area",
  data    = dat,
  backend = "TMB",
  mesh    = mesh,
  time    = "year"
)
```

## Testing

Tests use [`testthat`](https://testthat.r-lib.org/) and run automatically during `R CMD check`.

```r
devtools::test()
```
