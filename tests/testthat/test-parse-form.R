library(testthat)
library(dplyr)

# Source the function to be tested
source(here::here("R", "pre_proc.R"))

test_that("parse_form correctly parses formulas", {
  # 1. Create a sample data frame
  sample_data <- data.frame(
    bin = factor(c(1, 2, 1, 2)),
    gear = factor(c("A", "B", "A", "B")),
    area = factor(c("X", "X", "Y", "Y")),
    yy = c(2020, 2020, 2021, 2021),
    cvar1 = rnorm(4),
    tot_by_bin = rpois(4, 10),
    n = rpois(4, 100)
  )

  # 2. Test with brms backend
  parsed_brms <- parse_form(
    data = sample_data,
    backend = "brms",
    form = "gear + area + area:yy"
  )
  expect_s3_class(parsed_brms$form, "formula")
  expect_equal(
    as.character(parsed_brms$form),
    c("~", "tot_by_bin", "offset(log(n)) + 0 + bin + (1 | bin:gear) + (1 | bin:area) + (1 | bin:area:yy)")
  )

  # 3. Test with TMB backend
  parsed_tmb <- parse_form(
    data = sample_data,
    backend = "TMB",
    form = "gear + area"
  )
  expect_s3_class(parsed_tmb$form, "formula")
  expect_true(all(c("bin:gear", "bin:area") %in% names(parsed_tmb$data)))
})
