# ComPoM 0.4.0

## Changes

* New `step_comps()` and `plot_step_comps()`: refit with terms added one at a time and track a
  statistic of the annual composition (mean bin by default; `stat = 0.9` for the 90th percentile
  or a function; see `comp_quantile()`) for the catch and index components. TMB backend only.
* New `influ_comps()` and `plot_cdi_comps()`: the effect of each level on every bin, the
  distribution of the data over levels by year, and the influence of each term on the annual
  composition (log ratio by bin). `plot_cdi_comps()` follows Bentley et al. (2012): levels are
  shared between the coefficient and distribution panels, and years between the distribution and
  influence panels. `scale_comps(mode = "catch")` accepts `selectivity` to set
  terms to a reference.

* `scale_comps()` gains `mode = "index"` for the composition of the population that a CPUE index
  describes. Selectivity terms (`selectivity`, `reference`) and the `re_form` random effects
  are set to references, all other terms are predicted on a user-supplied `grid`, and rows are
  weighted by `weight` (equal if `NULL`). `remove_nuisance` is deprecated.
  `fit_model()` stores the fixed-effect variables in `ffx_vars`.

* `fit_model()` gains `spatiotemporal` (default `TRUE`). `spatiotemporal = FALSE` fits spatial
  fields only when a `mesh` is supplied.
* `cvar_smooth_plot()` and `spatial_field_plot()` no longer add a plot title.
* New `theme_compom()` (cowplot theme with faded major gridlines), used by all plotting
  functions. `spatial_field_plot()` uses the cividis scale and the smooth, ridge and composition
  plots use skyblue and orange.
* `scaled_comp_plot()` gains `comp_label` and `comp2_label` and shows a legend for the
  reference points. Its x axis has whole-number breaks.
* `spl = "2D"` with `backend = "TMB"` gave a non-positive-definite Hessian in all simulations
  tried (the bin margin of the tensor smooth is collinear with the bin intercepts). This is
  documented in `?parse_form`. Use `spl = "by_factor"` with `backend = "TMB"`.

# ComPoM 0.3.0

## Fixes

* `fit_model(backend = "brms")` failed with the default `dist`, because the
  default was the string `"poisson()"`. The default is now `"poisson"`.
* `data_prep()` no longer triggers the tidyselect external-vector deprecation
  warning.
* `post_pred_group()` works with `grp = NULL` and gives an informative error when the fit
  has no predictive draws (`add_preds = FALSE`).
* Package load no longer prints import-conflict warnings. Dependencies are
  imported selectively, and `recipes`, `stringr` and `posterior` are no longer
  dependencies.

## Changes

* `fit_model()` defaults `cores` to `min(chains, parallel::detectCores())` and `threads` to
  `NULL` (previously `32` and `16`).
* `fit_model()` warns when an sdmTMB fit has a Hessian that is not positive definite, and
  reports brms convergence with `message()` instead of `cat()`.
* `parse_form()` reports the formula with `message()` instead of printing it.
* `Fx_plot()`, `scaled_comp_plot()` and `scaled_ridge_plot()` gain an `xlab` argument. The
  default is the previous label, `"Length bin (cm)"`.
* `Selectivity()` is no longer part of the package. The simulation vignettes define it in
  `vignettes/vignette_helpers.R`.

## Additions

* Simulated example data `ex_comp` and `ex_catch`, used in all examples.
* Examples for all exported functions, and a documented return value for each.
* Tests of the covariate smooths, spatial field plots, `fit_model()`, `data_prep()` options,
  recovery of a known composition and the end-to-end workflow.
* MIT licence file, `NEWS.md`, `CONTRIBUTING.md`, citation, and a pkgdown site.
* GitHub Actions for `R CMD check` on several platforms, test coverage, vignette builds, the
  pkgdown site and releases.

# ComPoM 0.2.0

* Initial documented version with brms and sdmTMB backends, covariate smooths, spatial and
  spatio-temporal fields, and catch-weighted and CPUE-based scaling in `scale_comps()`.
