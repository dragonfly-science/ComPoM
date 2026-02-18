library(testthat)
library(dplyr)

# Source the function to be tested
source(here::here("R", "pre_proc.R"))

test_that("data_prep correctly processes compositional data", {
  # 1. Create a sample compositional data frame
  comp_data <- data.frame(
    year = rep(2020, 4),
    area = rep("A", 4),
    bin = c(10, 20, 10, 30),
    count = c(5, 10, 2, 8)
  )

  # 2. Run data_prep() on the sample data
  prepared_data <- data_prep(
    comp = comp_data,
    vars_for_grouping = c("year", "area"),
    bin_lab = "bin",
    count_var = "count"
  )

  # 3. Check that the output has the expected structure and values
  expect_true(all(c("year", "area", "bin", "tot_by_bin", "n") %in% names(prepared_data)))
  expect_equal(nrow(prepared_data), 3)
  expect_equal(prepared_data$tot_by_bin, c(7, 10, 8))
  expect_equal(prepared_data$n, c(25, 25, 25))
  expect_equal(prepared_data$bin, c(10, 20, 30))
})
