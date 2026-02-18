## Tests for pre-processing functions (data_prep, parse_form)
## Fixture objects (sim_comp, sim_dat) are created in helper-sim.R

# ── data_prep ─────────────────────────────────────────────────────────────────

test_that("data_prep returns expected columns", {
  expect_true(all(c("year", "area", "bin", "tot_by_bin", "n") %in% names(sim_dat)))
})

test_that("data_prep aggregates counts correctly", {
  # All bins present for each year x area combination
  expect_equal(nrow(sim_dat), 4L * 2L * 2L)  # 4 bins x 2 years x 2 areas
})

test_that("data_prep computes n as total counts per group", {
  totals <- sim_dat |>
    dplyr::group_by(year, area) |>
    dplyr::summarise(check = dplyr::first(n) == sum(tot_by_bin), .groups = "drop")
  expect_true(all(totals$check))
})

test_that("data_prep drops groups with zero total", {
  zero_comp <- sim_comp |> dplyr::mutate(count = 0L)
  result <- suppressWarnings(data_prep(zero_comp, c("year", "area"), "bin", "count"))
  expect_equal(nrow(result), 0L)
})

# ── parse_form ────────────────────────────────────────────────────────────────

test_that("parse_form returns a list with 'form' and 'data'", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", form = "year + area")
  )
  expect_type(pf, "list")
  expect_true(all(c("form", "data") %in% names(pf)))
})

test_that("parse_form formula is a formula object", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", form = "year + area")
  )
  expect_s3_class(pf$form, "formula")
})

test_that("parse_form formula contains offset(log(n)) for brms backend", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", form = "year + area")
  )
  # deparse() may return a multi-element character vector for long formulas
  expect_true(any(grepl("offset", deparse(pf$form))))
})

test_that("parse_form converts grouping columns to factors", {
  pf <- suppressMessages(
    parse_form(data = sim_dat, backend = "brms", form = "year + area")
  )
  expect_true(is.factor(pf$data$bin))
})
