library(testthat)
library(dplyr)

test_that("parse_form correctly parses formulas", {
  sample_data <- data.frame(
    bin        = factor(c(1, 2, 1, 2)),
    gear       = factor(c("A", "B", "A", "B")),
    area       = factor(c("X", "X", "Y", "Y")),
    yy         = factor(c(2020, 2020, 2021, 2021)),
    tot_by_bin = rpois(4, 10),
    n          = rpois(4, 100)
  )

  # brms: ffx_form = primary fixed effects, re_form = nuisance RE
  parsed_brms <- parse_form(
    data     = sample_data,
    backend  = "brms",
    ffx_form = "gear + area",
    re_form  = "yy"
  )
  expect_s3_class(parsed_brms$form, "formula")
  # formula must contain factor(gear):bin and factor(area):bin fixed effects
  form_str <- paste(deparse(parsed_brms$form), collapse = " ")
  expect_true(grepl("factor\\(gear\\):bin", form_str))
  expect_true(grepl("factor\\(area\\):bin", form_str))
  # and the nuisance (1|bin:yy) random effect
  expect_true(grepl("bin:yy", form_str))
  # re_vars should store the nuisance variable
  expect_equal(parsed_brms$re_vars, "yy")

  # TMB: interaction columns for nuisance RE should be created in data
  parsed_tmb <- parse_form(
    data     = sample_data,
    backend  = "TMB",
    ffx_form = "gear",
    re_form  = "area"
  )
  expect_s3_class(parsed_tmb$form, "formula")
  expect_true("bin:area" %in% names(parsed_tmb$data))
  expect_equal(parsed_tmb$re_vars, "area")
})
