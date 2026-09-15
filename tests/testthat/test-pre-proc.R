## Tests for pre-processing functions (data_prep, parse_form)
## Fixture objects (sim_comp, sim_dat) are created in helper-sim.R

# ── data_prep ─────────────────────────────────────────────────────────────────

test_that("data_prep returns expected columns", {
  expect_true(all(c("year", "area", "season", "bin", "tot_by_bin", "n") %in% names(sim_dat)))
})

test_that("data_prep aggregates counts correctly", {
  # 4 bins x 2 years x 2 areas x 2 seasons
  expect_equal(nrow(sim_dat), 4L * 2L * 2L * 2L)
})

test_that("data_prep computes n as total counts per group", {
  totals <- sim_dat |>
    dplyr::group_by(year, area, season) |>
    dplyr::summarise(check = dplyr::first(n) == sum(tot_by_bin), .groups = "drop")
  expect_true(all(totals$check))
})

test_that("data_prep drops groups with zero total", {
  zero_comp <- sim_comp |> dplyr::mutate(count = 0L)
  result <- suppressWarnings(data_prep(zero_comp, c("year", "area", "season"), "bin", "count"))
  expect_equal(nrow(result), 0L)
})

# ── parse_form ────────────────────────────────────────────────────────────────

test_that("parse_form returns a list with 'form', 'data', 're_vars'", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_type(pf, "list")
  expect_true(all(c("form", "data", "re_vars") %in% names(pf)))
})

test_that("parse_form formula is a formula object", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_s3_class(pf$form, "formula")
})

test_that("parse_form formula contains offset(log(n)) for brms backend", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_true(any(grepl("offset", deparse(pf$form))))
})

test_that("parse_form formula contains factor(year):bin fixed effect", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_true(any(grepl("factor\\(year\\)", deparse(pf$form))))
})

test_that("parse_form formula contains (1|bin:season) nuisance RE", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_true(any(grepl("bin:season", deparse(pf$form))))
})

test_that("parse_form re_vars stores nuisance variable names", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_equal(pf$re_vars, "season")
})

test_that("parse_form re_vars is NULL when re_form is NULL", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = NULL)
  )
  expect_null(pf$re_vars)
})

test_that("parse_form converts grouping columns to factors", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", ffx_form = "year", re_form = "season")
  )
  expect_true(is.factor(pf$data$bin))
})

# ── fit_model stores re_vars ──────────────────────────────────────────────────

test_that("fit_model stores re_vars in ComPoM object [brms]", {
  expect_equal(sim_fit_brms$re_vars, "season")
})

test_that("fit_model stores re_vars in ComPoM object [TMB non-spatial]", {
  expect_equal(sim_fit_tmb$re_vars, "season")
})

test_that("fit_model re_vars is NULL when re_form not supplied [brms]", {
  fit_no_re <- suppressMessages(suppressWarnings(
    fit_model(
      ffx_form     = "year",
      re_form      = NULL,
      data         = sim_dat,
      dist         = poisson(),
      backend      = "brms",
      brms_backend = "rstan",
      chains       = 1, cores = 1, threads = 1,
      iter = 400, warmup = 200, thin = 1, refresh = 0
    )
  ))
  expect_null(fit_no_re$re_vars)
})
