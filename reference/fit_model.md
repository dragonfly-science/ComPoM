# Fit the Poisson-multinomial model

Builds the model formula with
[`parse_form()`](https://dragonfly-science.github.io/ComPoM/reference/parse_form.md)
and fits it with either `brms` (Bayesian, `backend = "brms"`) or
`sdmTMB` (maximum likelihood, `backend = "TMB"`). Counts per stratum and
bin are modelled as Poisson with `log(n)` as offset, so that the implied
within-stratum proportions follow a multinomial-logit model (see
[`vignette("ComPoM")`](https://dragonfly-science.github.io/ComPoM/articles/ComPoM.md)
for the model equations).

## Usage

``` r
fit_model(
  ffx_form,
  re_form = NULL,
  cvars = NULL,
  spl = "2D",
  knots = "c(5,10)",
  data = NULL,
  dist = "poisson",
  chains = 4,
  threads = NULL,
  cores = min(chains, parallel::detectCores()),
  backend = "brms",
  brms_backend = "rstan",
  iter = 1000,
  warmup = 500,
  mtd = 10,
  ad = 0.8,
  thin = 2,
  refresh = 10,
  add_preds = FALSE,
  mesh = NULL,
  time = NULL,
  share_spatial_sd = TRUE,
  share_spatiotemporal_sd = TRUE,
  spatiotemporal = TRUE,
  coords = c("x", "y"),
  priors = NULL,
  newton_loops = 1
)
```

## Arguments

- ffx_form:

  A character string of primary fixed-effect variables (RHS), e.g.
  `"Year"`. Passed to
  [`parse_form()`](https://dragonfly-science.github.io/ComPoM/reference/parse_form.md).

- re_form:

  A character string of nuisance random-effect variables, e.g.
  `"Season + Gear"`. These become `(1|bin:var)` IID intercepts. `NULL` =
  none.

- cvars:

  Continuous covariate(s) to be fitted as smooths with bin; see `spl`.

- spl:

  Smooth type for `cvars`: `"2D"`, `"by_factor"` or `"by"`; see
  [`parse_form()`](https://dragonfly-science.github.io/ComPoM/reference/parse_form.md).
  `"2D"` is unreliable with `backend = "TMB"` (non-positive-definite
  Hessian); use `"by_factor"` there.

- knots:

  Basis dimension(s) for the `cvars` smooth, as a character string (e.g.
  `"4"` or `"c(5,10)"` for `spl = "2D"`).

- data:

  Prepared data frame from
  [`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md).

- dist:

  Response family: a family name (default `"poisson"`) or a family
  object such as [`poisson()`](https://rdrr.io/r/stats/family.html).
  With `backend = "brms"` it is passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) as
  `family`. With `backend = "TMB"` it must be `poisson`, `nbinom2` or
  `nbinom1` (sdmTMB). The negative binomial families allow variance
  between bins in excess of the Poisson variance. The offset is kept,
  but the equivalence to a multinomial likelihood then holds only
  approximately (see
  [`vignette("stratum_intercept")`](https://dragonfly-science.github.io/ComPoM/articles/stratum_intercept.md)).

- chains, cores, threads:

  Number of MCMC chains, cores and within-chain threads (brms backend
  only). By default, one core is used per chain (up to the number of
  available cores) and no within-chain threading.

- backend:

  the backend used for model fitting: `"brms"` or `"TMB"`.

- brms_backend:

  Stan interface used by brms (`"rstan"` or `"cmdstanr"`).

- iter, warmup, thin:

  MCMC iterations, warm-up iterations and thinning (brms). For the TMB
  backend these only set the number of simulation draws used downstream,
  `nsim = (iter - warmup) / thin`.

- mtd, ad:

  Stan `max_treedepth` and `adapt_delta` (brms only).

- refresh:

  Stan progress refresh interval (brms only).

- add_preds:

  Add posterior predictive draws for the fitted data, required by
  [`post_pred_group()`](https://dragonfly-science.github.io/ComPoM/reference/post_pred_group.md).
  Can be slow.

- mesh:

  An sdmTMB mesh object created with make_mesh() (TMB backend only)

- time:

  Column name for time (year) variable (TMB backend only)

- share_spatial_sd:

  Share spatial SD across bins (TMB backend only)

- share_spatiotemporal_sd:

  Share spatiotemporal SD across bins (TMB backend only)

- spatiotemporal:

  Include spatio-temporal fields (one per bin and time step) in addition
  to the spatial fields (TMB backend with a `mesh` only). `FALSE` fits
  spatial fields only; `time` is still required, as it is used for
  prediction.

- coords:

  Character vector of coordinate column names (TMB backend only)

- priors:

  An
  [`sdmTMBpriors()`](https://sdmTMB.github.io/sdmTMB/reference/priors.html)
  object for regularisation (TMB non-spatial backend only; currently not
  passed to the spatial fit). Useful when the design matrix is sparse
  (e.g. `normal(0, 2)` on fixed effects).

- newton_loops:

  Number of extra Newton optimisation steps after initial fit (TMB
  non-spatial only). Default 1. Increase to 2-3 for difficult
  convergence. For spatial fits use
  [`sdmTMB::run_extra_optimization()`](https://sdmTMB.github.io/sdmTMB/reference/run_extra_optimization.html)
  on `$model`.

## Value

An object of class `ComPoM`: a list with the fitted model (`model`, a
`brmsfit` or `sdmTMB` object), the model data (`data`), the number of
draws used downstream (`nsim`), posterior predictive draws (`preds`, or
`NULL`), and the settings needed by the post-processing functions
(`mesh`, `spatial_varying`, `time`, `coords`, `ffx_vars`, `re_vars`,
`cvars`).

## Details

With `backend = "TMB"` and a `mesh`, the model adds one spatial Gaussian
random field per bin and one spatio-temporal field per bin and time step
(set up with
[`sdmTMB::make_category_svc()`](https://sdmTMB.github.io/sdmTMB/reference/make_category_svc.html)).
Set `spatiotemporal = FALSE` to fit the spatial fields only.

With `backend = "brms"`, fixed effects receive a `student_t(10, 0, 1)`
prior. With `backend = "TMB"`, uncertainty in downstream functions is
obtained from `nsim` draws from the approximate joint sampling
distribution of the parameters.

## See also

[`data_prep()`](https://dragonfly-science.github.io/ComPoM/reference/data_prep.md),
[`scale_comps()`](https://dragonfly-science.github.io/ComPoM/reference/scale_comps.md),
[`Fx_plot()`](https://dragonfly-science.github.io/ComPoM/reference/Fx_plot.md),
[`vignette("ComPoM")`](https://dragonfly-science.github.io/ComPoM/articles/ComPoM.md)

## Examples

``` r
dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")

# maximum likelihood fit: Year as effect of interest, Season and Vessel as
# nuisance random effects
fit <- fit_model(
  ffx_form = "Year", re_form = "Season + Vessel",
  data = dat, backend = "TMB", iter = 300, warmup = 100, thin = 2
)
#> Model formula: tot_by_bin ~ 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)
class(fit)
#> [1] "ComPoM"

# \donttest{
# Bayesian fit with brms (compiles a Stan model; slow)
fit_b <- fit_model(
  ffx_form = "Year", re_form = "Season + Vessel", data = dat,
  backend = "brms", chains = 2, iter = 400, warmup = 200, thin = 1
)
#> Model formula: tot_by_bin ~ offset(log(n)) + 0 + bin + factor(Year):bin + (1 | bin:Season) + (1 | bin:Vessel)
#> Compiling Stan program...
#> Start sampling
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Converged (all Rhat < 1.05): TRUE
# }
```
