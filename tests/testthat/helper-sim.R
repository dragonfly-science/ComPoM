## ── Shared test fixtures ──────────────────────────────────────────────────────
## testthat automatically sources any file named helper-*.R before running tests.

library(ComPoM)
library(dplyr)
library(tidyr)
library(sdmTMB)

# ── 1. Minimal simulated compositional dataset ────────────────────────────────
set.seed(42)

sim_comp <- expand.grid(
  year   = factor(2020:2021),
  area   = factor(c("A", "B")),
  season = factor(c("S1", "S2")),
  bin    = c(10, 20, 30, 40)
) |>
  mutate(count = rpois(n(), lambda = 10 + as.integer(bin) / 5))

# ── 2. Prepare data ───────────────────────────────────────────────────────────
sim_dat <- data_prep(
  comp              = sim_comp,
  vars_for_grouping = c("year", "area", "season"),
  bin_lab           = "bin",
  count_var         = "count"
)

# ── 3a. Fit a minimal brms model (1 chain, 200 posterior draws) ───────────────
#        ffx_form = primary effect; re_form = nuisance RE to be dropped at scale
# The brms fit needs rstan (Suggests) and a C++ toolchain. If it is unavailable, the
# fixture is NULL and brms tests are skipped with skip_if_no_brms(). The fit is also
# skipped when the Stan model cannot be compiled (e.g. rstan's sink() error under testthat on
# Windows), so that a failure in the compiler does not abort the whole test run.
sim_fit_brms <- if (!requireNamespace("rstan", quietly = TRUE)) NULL else tryCatch(suppressMessages(
  suppressWarnings(
    fit_model(
      ffx_form     = "year",
      re_form      = "season",
      data         = sim_dat,
      dist         = poisson(),
      backend      = "brms",
      brms_backend = "rstan",
      chains       = 1,
      cores        = 1,
      threads      = 1,
      iter         = 400,
      warmup       = 200,
      thin         = 1,
      refresh      = 0,
      add_preds    = TRUE
    )
  )
), error = function(e) {
  message("brms fixture not available: ", conditionMessage(e))
  NULL
})

# keep old alias so existing tests still pass
sim_fit <- sim_fit_brms

# ── 3b. Fit a minimal sdmTMB non-spatial model ────────────────────────────────
sim_fit_tmb <- suppressMessages(
  suppressWarnings(
    fit_model(
      ffx_form  = "year",
      re_form   = "season",
      data      = sim_dat,
      backend   = "TMB",
      iter      = 400,
      warmup    = 200,
      thin      = 1,
      add_preds = TRUE
    )
  )
)

# ── 3c. Fit a minimal sdmTMB spatial model ────────────────────────────────────
# Strata are sampled at distinct random locations in every year. With one
# bin-specific spatio-temporal field per bin and year, a design with few repeated
# locations (e.g. 4 locations x 2 years) is not identifiable and the Hessian is
# not positive definite, which makes simulate()/predict(nsim) fail.
set.seed(1)
sim_strata_spatial <- expand.grid(year = factor(2020:2022), k = 1:10) |>
  mutate(
    x      = runif(n(), 0, 100),
    y      = runif(n(), 0, 100),
    area   = factor(ifelse(x < 50, "A", "B")),
    season = factor(sample(c("S1", "S2"), n(), replace = TRUE)),
    loc    = factor(row_number())
  )

# Spatial gradient, year trend and season effect that differ by bin. Counts with no such
# signal give boundary estimates (field SDs -> 0 or Inf) and a non-positive-definite
# joint precision matrix.
sim_comp_spatial <- merge(sim_strata_spatial, data.frame(bin = c(10, 20, 30, 40))) |>
  mutate(
    b     = as.integer(factor(bin)),
    eta   = log(40) + 0.02 * (x - 50) * (b - 2.5) + 0.02 * (y - 50) * (2.5 - b) / 2 +
            0.3 * (as.integer(year) - 2) * (b - 2.5) / 2 +
            0.4 * (season == "S2") * (b - 2.5) / 1.5,
    count = rpois(n(), exp(eta))
  )

# x and y are factors after data_prep(); convert back to numeric for the mesh
sim_dat_spatial <- data_prep(
  comp              = sim_comp_spatial,
  vars_for_grouping = c("year", "area", "season", "loc", "x", "y"),
  bin_lab           = "bin",
  count_var         = "count"
) |>
  mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))

sim_mesh <- suppressMessages(
  sdmTMB::make_mesh(sim_dat_spatial, xy_cols = c("x", "y"), cutoff = 40)
)

sim_fit_tmb_spatial <- suppressMessages(
  suppressWarnings(
    fit_model(
      ffx_form  = "year",
      re_form   = "season",
      data      = sim_dat_spatial,
      backend   = "TMB",
      iter      = 400,
      warmup    = 200,
      thin      = 1,
      mesh      = sim_mesh,
      time      = "year",
      coords    = c("x", "y"),
      add_preds = TRUE
    )
  )
)

# ── 4. Minimal scale_df for scale_comps() ────────────────────────────────────
# season is the nuisance RE in the fixtures; catch weighting
# (remove_nuisance = FALSE) requires it in scale_df and grps
sim_scale_df <- expand.grid(
  year   = factor(2020:2021),
  area   = factor(c("A", "B")),
  season = factor(c("S1", "S2"))
) |>
  mutate(catch = rep(c(100, 150, 200, 250), 2))

sim_scale_df_spatial <- sim_strata_spatial |>
  dplyr::select(year, area, season, x, y) |>
  mutate(catch = rep(c(100, 150, 200, 250, 120), length.out = n()))

skip_if_no_brms <- function() {
  testthat::skip_if(is.null(sim_fit_brms), "rstan is not installed or the Stan model did not compile")
}

# two fixed-effect terms (year, area) and a nuisance random effect (season)
sim_fit_two <- suppressMessages(suppressWarnings(
  fit_model(ffx_form = "year + area", re_form = "season", data = sim_dat, backend = "TMB",
            iter = 300, warmup = 100, thin = 1)
))
