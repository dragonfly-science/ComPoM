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
sim_fit_brms <- suppressMessages(
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
)

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
group_coords <- tibble::tribble(
  ~year,        ~area, ~season, ~x,  ~y,
  factor(2020), "A",   "S1",    0,   0,
  factor(2020), "A",   "S2",    0,   0,
  factor(2020), "B",   "S1",   10,   0,
  factor(2020), "B",   "S2",   10,   0,
  factor(2021), "A",   "S1",    0,  10,
  factor(2021), "A",   "S2",    0,  10,
  factor(2021), "B",   "S1",   10,  10,
  factor(2021), "B",   "S2",   10,  10
) |>
  mutate(area = factor(area), season = factor(season))

sim_dat_spatial <- sim_dat |>
  left_join(group_coords, by = c("year", "area", "season"))

sim_mesh <- suppressMessages(
  make_mesh(sim_dat_spatial, xy_cols = c("x", "y"), cutoff = 3)
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
sim_scale_df <- expand.grid(
  year  = factor(2020:2021),
  area  = factor(c("A", "B"))
) |>
  mutate(catch = c(100, 150, 200, 250))

sim_scale_df_spatial <- sim_scale_df |>
  mutate(x = c(0, 0, 10, 10), y = c(0, 10, 0, 10))
