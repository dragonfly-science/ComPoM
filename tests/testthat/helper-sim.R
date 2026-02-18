## ── Shared test fixtures ──────────────────────────────────────────────────────
## testthat automatically sources any file named helper-*.R before running tests.

library(ComPoM)
library(dplyr)
library(tidyr)

# ── 1. Minimal simulated compositional dataset ────────────────────────────────
set.seed(42)

sim_comp <- expand.grid(
  year  = factor(2020:2021),
  area  = factor(c("A", "B")),
  bin   = c(10, 20, 30, 40)
) |>
  mutate(count = rpois(n(), lambda = 10 + as.integer(bin) / 5))

# ── 2. Prepare data ───────────────────────────────────────────────────────────
sim_dat <- data_prep(
  comp              = sim_comp,
  vars_for_grouping = c("year", "area"),
  bin_lab           = "bin",
  count_var         = "count"
)

# ── 3. Fit a minimal brms model (1 chain, 200 posterior draws) ────────────────
sim_fit <- suppressMessages(
  suppressWarnings(
    fit_model(
      form         = "year + area",
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

# ── 4. Minimal scale_df for scale_comps() ────────────────────────────────────
sim_scale_df <- expand.grid(
  year  = factor(2020:2021),
  area  = factor(c("A", "B"))
) |>
  mutate(catch = c(100, 150, 200, 250))
