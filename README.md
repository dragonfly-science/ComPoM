# ComPoM

<!-- badges: start -->
[![R-CMD-check](https://github.com/dragonfly-science/ComPoM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dragonfly-science/ComPoM/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/dragonfly-science/ComPoM/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/dragonfly-science/ComPoM/actions/workflows/test-coverage.yaml)
[![Codecov test coverage](https://codecov.io/gh/dragonfly-science/ComPoM/graph/badge.svg)](https://app.codecov.io/gh/dragonfly-science/ComPoM)
[![pkgdown](https://github.com/dragonfly-science/ComPoM/actions/workflows/pkgdown.yaml/badge.svg)](https://dragonfly-science.github.io/ComPoM/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
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

Documentation, function reference and vignettes are at
<https://dragonfly-science.github.io/ComPoM/>.

### Requirements

- R >= 4.1.
- The `sdmTMB` backend (`backend = "TMB"`) needs [sdmTMB](https://sdmtmb.github.io/sdmTMB/),
  which compiles C++ code on installation and needs a C++ toolchain (Rtools on Windows, Xcode
  command line tools on macOS).
- The `brms` backend (`backend = "brms"`) fits Stan models and needs a Stan interface and a C++
  toolchain. Install `rstan` (the default, `brms_backend = "rstan"`) or
  [cmdstanr](https://mc-stan.org/cmdstanr/) (`brms_backend = "cmdstanr"`). See the
  [brms](https://paul-buerkner.github.io/brms/) and Stan installation notes.
- The vignettes need the suggested packages (`tidyverse`, `patchwork`, `sf`, `terra`, `fmesher`,
  `MASS`).

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

The package includes simulated sampling data (`ex_comp`, one row per stratum and length bin) and
the catch of the strata (`ex_catch`), see `?ex_data`.

```r
library(ComPoM)

# 1. Prepare data: one row per stratum x bin. Include every variable that is
#    used later (fixed effects, nuisance effects, covariates, coordinates).
dat <- data_prep(
  comp              = ex_comp,
  vars_for_grouping = c("Year", "Season", "Vessel", "trip"),
  bin_lab           = "length",
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
#     ex_catch has one row per stratum with its catch; all nuisance variables
#     are predicted at each stratum's level.
scaled <- scale_comps(
  scale_df = ex_catch,
  predvar  = "catch",
  fit      = fit,
  grps     = c("Year", "Season", "Vessel")
)
scaled_comp_plot(scaled, grps = "Year")

# 4b. CPUE-based scaling: composition representative of the population,
#     standardised for the nuisance variables (removed from the prediction).
#     One weight (e.g. a standardised CPUE) per year.
cpue <- data.frame(Year = factor(2015:2020), cpue = 1)
scaled_pop <- scale_comps(
  scale_df        = cpue,
  predvar         = "cpue",
  fit             = fit,
  grps            = "Year",
  remove_nuisance = TRUE
)
scaled_ridge_plot(scaled_pop, grps = "Year", plot_y = "Year")
```

With `backend = "brms"`, `fit_model()` samples with Stan (`chains`, `iter`, `warmup`, `thin`
control the sampling), and the draws are posterior draws. With `backend = "TMB"` the draws are
simulated from the approximate joint distribution of the estimates (`nsim = (iter - warmup) /
thin`), so `iter`, `warmup` and `thin` only set the number of draws. Check the
`sdmTMB` fit with `sdmTMB::sanity(fit$model)`, and see the warning raised by `fit_model()` when
the Hessian is not positive definite.

### Continuous covariate smooth

`sst` is a variable in `ex_comp`. It has to be included in `vars_for_grouping`, and be numeric.
In `scale_comps()`, covariates are averaged within each stratum of `scale_df`, so set them to a
reference value in `scale_df` to standardise to that value.

```r
dat_sst <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip", "sst"), "length", "count")
dat_sst$sst <- as.numeric(as.character(dat_sst$sst))
fit_sst <- fit_model(
  ffx_form = "Year", re_form = "Season", data = dat_sst, backend = "TMB",
  cvars = "sst", spl = "by_factor", knots = "4"
)
cvar_smooth_plot(fit_sst, "sst")
```

### Spatio-temporal model (sdmTMB backend)

```r
library(sdmTMB)

dat_sp <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip", "x", "y"), "length", "count")
# coordinates are factors after data_prep(); convert to numeric
dat_sp <- dplyr::mutate(dat_sp, x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
mesh <- make_mesh(dat_sp, xy_cols = c("x", "y"), cutoff = 25)

fit_sp <- fit_model(
  ffx_form = "Year",
  re_form  = "Season + Vessel",
  data     = dat_sp,
  backend  = "TMB",
  mesh     = mesh,
  time     = "Year",
  coords   = c("x", "y")
)
spatial_field_plot(fit_sp)
```

For spatial models, the scaling data passed to `scale_comps()` must contain the coordinates, and
the time column must have the same class (e.g. factor) as in the model data. Spatial and
spatio-temporal fields are estimated for every bin, so they need data at many locations per
year. With few locations the model is not identifiable and the Hessian is not positive definite.

## Vignettes

- `ComPoM`: model equations (basic and spatio-temporal), workflow and a simulated example.
- `vessel_bias_simulation`: correction for vessel-level sampling bias (brms and sdmTMB).
- `spatial_bias_simulation`: correction for spatially biased sampling, and covariate smooths.

The vignettes are slow, because they fit models (the spatial vignette takes about 10 minutes).
They are not built by `R CMD check` in CI. The `vignettes` workflow renders them, and the
rendered versions are on the [package website](https://dragonfly-science.github.io/ComPoM/).

## Testing

Tests use [testthat](https://testthat.r-lib.org/):

```r
devtools::test()
```

The tests include recovery of a known composition from simulated data and an end-to-end run of the
workflow above. The brms tests are skipped if `rstan` is not installed.

GitHub Actions run `R CMD check` (Linux, macOS, Windows), test coverage with
[covr](https://covr.r-lib.org/), the vignettes, and the pkgdown site on every push and pull
request. Coverage is reported on [Codecov](https://app.codecov.io/gh/dragonfly-science/ComPoM).
To compute coverage locally, run

```r
covr::report(covr::package_coverage())
```

## Contributing and citation

See [CONTRIBUTING.md](CONTRIBUTING.md) for bug reports and development. Cite the package with
`citation("ComPoM")`. Release notes are in [NEWS.md](NEWS.md). ComPoM is released under the
[MIT licence](LICENSE.md).
