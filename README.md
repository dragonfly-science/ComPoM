# ComPoM

<!-- badges: start -->
[![R-CMD-check](https://github.com/dragonfly-science/ComPoM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dragonfly-science/ComPoM/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/dragonfly-science/ComPoM/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/dragonfly-science/ComPoM/actions/workflows/test-coverage.yaml)
[![Codecov test coverage](https://codecov.io/gh/dragonfly-science/ComPoM/graph/badge.svg)](https://app.codecov.io/gh/dragonfly-science/ComPoM)
<!-- badges: end -->

ComPoM fits Poisson-multinomial models to compositional count data (e.g. length or age
frequencies from observer sampling) and uses them to standardise compositions for
non-representative sampling and scale them to total catch.

Counts per bin and sampling stratum are modelled as Poisson with the stratum sample size as
offset. All covariate effects are bin-specific, so that the implied within-stratum proportions
follow a multinomial-logit model:

- effects of interest (e.g. year) enter as bin-specific fixed effects (`ffx_form`);
- nuisance sampling factors (e.g. vessel, season) enter as bin-specific IID random effects
  (`re_form`);
- continuous covariates (e.g. sea-surface temperature) enter as smooths over covariate and bin
  (`cvars`, `spl`);
- with the sdmTMB backend and a mesh, bin-specific spatial and spatio-temporal Gaussian random
  fields are added.

Models are fitted with [brms](https://paul-buerkner.github.io/brms/) (Bayesian) or
[sdmTMB](https://sdmtmb.github.io/sdmTMB/) (maximum likelihood). The model equations are given in
`vignette("ComPoM")`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("dragonfly-science/ComPoM")
```

Or from a local clone, in the package root:

```r
devtools::install()
```

## Workflow

| Step | Function | Purpose |
|---|---|---|
| 1 | `data_prep()` | Aggregate counts to strata x bins, add zero counts and stratum totals |
| 2 | `fit_model()` | Build the formula (via `parse_form()`) and fit with brms or sdmTMB |
| 3 | `post_pred_group()` | Posterior predictive check by group (needs `add_preds = TRUE`) |
| 3 | `Fx_plot()` | Bin-specific effects of a factor, as multipliers |
| 3 | `cvar_smooth_plot()` | Smooth effect of a continuous covariate by bin |
| 3 | `spatial_field_plot()` | Spatial field by bin (spatial sdmTMB fits) |
| 4 | `scale_comps()` | Predict numbers-at-bin draws for strata with known catch or CPUE |
| 4 | `scaled_comp_plot()`, `scaled_ridge_plot()` | Aggregate draws to compositions and plot |

Models are usually developed in increasing complexity: fixed effects only, then nuisance random
effects, then covariate smooths, then spatial fields.

## Basic usage

```r
library(ComPoM)

# 1. Prepare data: one row per stratum x bin. Include every variable that is
#    used later (fixed effects, nuisance effects, covariates, coordinates).
dat <- data_prep(
  comp              = my_data,
  vars_for_grouping = c("Year", "Season", "Vessel"),
  bin_lab           = "length_bin",
  count_var         = "count"
)

# 2. Fit: Year as effect of interest, Season and Vessel as nuisance effects
fit <- fit_model(
  ffx_form  = "Year",
  re_form   = "Season + Vessel",
  data      = dat,
  backend   = "TMB",     # or "brms"
  add_preds = TRUE
)

# 3. Diagnostics and effects
post_pred_group(fit, grp = "Year")
Fx_plot(fit, grp = "Year")

# 4a. Catch weighting: composition of the total catch of the fleet.
#     catch_data has one row per stratum (Year x Season x Vessel) with its
#     catch; all nuisance variables are predicted at each stratum's level.
scaled <- scale_comps(
  scale_df = catch_data,
  predvar  = "catch",
  fit      = fit,
  grps     = c("Year", "Season", "Vessel")
)
scaled_comp_plot(scaled, grps = "Year")

# 4b. CPUE-based scaling: composition representative of the population,
#     standardised for the nuisance variables (removed from the prediction).
#     cpue_data has a predicted CPUE or abundance per stratum (e.g. Year).
scaled_pop <- scale_comps(
  scale_df        = cpue_data,
  predvar         = "cpue",
  fit             = fit,
  grps            = "Year",
  remove_nuisance = TRUE
)
```

### Continuous covariate smooth

```r
fit_sst <- fit_model(
  ffx_form = "Year", re_form = "Season", data = dat, backend = "TMB",
  cvars = "sst", spl = "by_factor", knots = "4"
)
cvar_smooth_plot(fit_sst, "sst")
```

### Spatio-temporal model (sdmTMB backend)

```r
library(sdmTMB)

# coordinates are factors after data_prep(); convert to numeric
dat <- dplyr::mutate(dat, x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
mesh <- make_mesh(dat, xy_cols = c("x", "y"), cutoff = 10)

fit_sp <- fit_model(
  ffx_form = "Year",
  re_form  = "Season",
  data     = dat,
  backend  = "TMB",
  mesh     = mesh,
  time     = "Year",
  coords   = c("x", "y")
)
spatial_field_plot(fit_sp)
```

For spatial models, the scaling data passed to `scale_comps()` must contain the coordinates, and
the time column must have the same class (e.g. factor) as in the model data.

## Vignettes

- `ComPoM`: model equations (basic and spatio-temporal), workflow and a simulated example.
- `vessel_bias_simulation`: correction for vessel-level sampling bias (brms and sdmTMB).
- `spatial_bias_simulation`: correction for spatially biased sampling, and covariate smooths.

## Testing

Tests use [testthat](https://testthat.r-lib.org/):

```r
devtools::test()
```

GitHub Actions run `R CMD check` (without vignettes) and test coverage with
[covr](https://covr.r-lib.org/) on every push and pull request. Coverage is reported on
[Codecov](https://app.codecov.io/gh/dragonfly-science/ComPoM). To compute coverage locally, run

```r
covr::report(covr::package_coverage())
```
