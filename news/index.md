# Changelog

## ComPoM 0.3.0

### Fixes

- `fit_model(backend = "brms")` failed with the default `dist`, because
  the default was the string `"poisson()"`. The default is now
  `"poisson"`.
- [`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md)
  no longer triggers the tidyselect external-vector deprecation warning.
- [`post_pred_group()`](https://dragonfly-science.github.io/ComPoM/reference/post_pred_group.md)
  works with `grp = NULL` and gives an informative error when the fit
  has no predictive draws (`add_preds = FALSE`).
- Package load no longer prints import-conflict warnings. Dependencies
  are imported selectively, and `recipes`, `stringr` and `posterior` are
  no longer dependencies.

### Changes

- [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md)
  gains `mode = "index"` for the composition of the population that a
  CPUE index describes. Catchability terms (`catchability`, `reference`)
  and the `re_form` random effects are set to references, all other
  terms are predicted on a user-supplied `grid`, and rows are weighted
  by `weight` (equal if `NULL`). `remove_nuisance` is deprecated.
  [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
  stores the fixed-effect variables in `ffx_vars`.

- [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
  gains `spatiotemporal` (default `TRUE`). `spatiotemporal = FALSE` fits
  spatial fields only when a `mesh` is supplied.

- [`cvar_smooth_plot()`](https://dragonfly-science.github.io/ComPoM/reference/cvar_smooth_plot.md)
  and
  [`spatial_field_plot()`](https://dragonfly-science.github.io/ComPoM/reference/spatial_field_plot.md)
  no longer add a plot title.

- New
  [`theme_compom()`](https://dragonfly-science.github.io/ComPoM/reference/theme_compom.md)
  (cowplot theme with faded major gridlines), used by all plotting
  functions.
  [`spatial_field_plot()`](https://dragonfly-science.github.io/ComPoM/reference/spatial_field_plot.md)
  uses the cividis scale and the smooth, ridge and composition plots use
  skyblue and orange.

- [`scaled_comp_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_comp_plot.md)
  gains `comp_label` and `comp2_label` and shows a legend for the
  reference points. Its x axis has whole-number breaks.

- `spl = "2D"` with `backend = "TMB"` gave a non-positive-definite
  Hessian in all simulations tried (the bin margin of the tensor smooth
  is collinear with the bin intercepts). This is documented in
  [`?parse_form`](https://dragonfly-science.github.io/ComPoM/reference/parse_form.md).
  Use `spl = "by_factor"` with `backend = "TMB"`.

- [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
  defaults `cores` to `min(chains, parallel::detectCores())` and
  `threads` to `NULL` (previously `32` and `16`).

- [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md)
  warns when an sdmTMB fit has a Hessian that is not positive definite,
  and reports brms convergence with
  [`message()`](https://rdrr.io/r/base/message.html) instead of
  [`cat()`](https://rdrr.io/r/base/cat.html).

- [`parse_form()`](https://dragonfly-science.github.io/ComPoM/reference/parse_form.md)
  reports the formula with
  [`message()`](https://rdrr.io/r/base/message.html) instead of printing
  it.

- [`Fx_plot()`](https://dragonfly-science.github.io/ComPoM/reference/Fx_plot.md),
  [`scaled_comp_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_comp_plot.md)
  and
  [`scaled_ridge_plot()`](https://dragonfly-science.github.io/ComPoM/reference/scaled_ridge_plot.md)
  gain an `xlab` argument. The default is the previous label,
  `"Length bin (cm)"`.

- `Selectivity()` is no longer part of the package. The simulation
  vignettes define it in `vignettes/vignette_helpers.R`.

### Additions

- Simulated example data `ex_comp` and `ex_catch`, used in all examples.
- Examples for all exported functions, and a documented return value for
  each.
- Tests of the covariate smooths, spatial field plots,
  [`fit_model()`](https://dragonfly-science.github.io/ComPoM/reference/fit_model.md),
  [`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md)
  options, recovery of a known composition and the end-to-end workflow.
- MIT licence file, `NEWS.md`, `CONTRIBUTING.md`, citation, and a
  pkgdown site.
- GitHub Actions for `R CMD check` on several platforms, test coverage,
  vignette builds, the pkgdown site and releases.

## ComPoM 0.2.0

- Initial documented version with brms and sdmTMB backends, covariate
  smooths, spatial and spatio-temporal fields, and catch-weighted and
  CPUE-based scaling in
  [`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md).
