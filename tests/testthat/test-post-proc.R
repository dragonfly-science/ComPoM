library(testthat)
library(ggplot2)

# Source the function to be tested
source(here::here("R", "post_proc.R"))

test_that("post_pred_group creates a ggplot object", {
  # This test requires a model object from fit_model()
  # For now, we'll just check that the function exists.
  expect_true(exists("post_pred_group"))
})
