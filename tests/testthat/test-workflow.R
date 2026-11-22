## End-to-end workflow on the bundled example data (the workflow shown in the README).
## sdmTMB backend only, so these tests are fast and need no Stan toolchain.

strata <- c("Year", "Season", "Vessel", "trip")

test_that("bundled example data have the documented structure", {
  expect_true(all(c("Year", "Season", "Vessel", "trip", "x", "y", "sst", "length", "count")
                  %in% names(ex_comp)))
  expect_equal(nrow(ex_comp), 144 * 6)
  expect_equal(nrow(ex_catch), 144)
  expect_true("catch" %in% names(ex_catch))
})

test_that("non-spatial workflow runs from data_prep to plots", {
  dat <- data_prep(ex_comp, strata, "length", "count")
  fit <- suppressMessages(fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                                    iter = 300, warmup = 100, thin = 2, add_preds = TRUE))
  expect_true(isTRUE(fit$model$sd_report$pdHess))
  expect_s3_class(post_pred_group(fit, grp = "Year"), "ggplot")
  expect_s3_class(Fx_plot(fit, grp = "Year"), "ggplot")

  sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
  expect_equal(length(unique(sc$.draw)), fit$nsim)
  expect_s3_class(scaled_comp_plot(sc, grps = "Year"), "ggplot")
  expect_s3_class(scaled_ridge_plot(sc, grps = "Year", plot_y = "Year"), "ggplot")

  # overlay of raw compositions
  raw <- dplyr::rename(ex_comp, bin = length)
  expect_s3_class(scaled_comp_plot(sc, grps = "Year", comp_are = raw, cvar = "count"), "ggplot")
})

test_that("post_pred_group pools strata when grp is NULL and errors without draws", {
  dat <- data_prep(ex_comp, strata, "length", "count")
  fit <- suppressMessages(fit_model("Year", "Season", data = dat, backend = "TMB",
                                    iter = 300, warmup = 100, thin = 2, add_preds = TRUE))
  expect_s3_class(post_pred_group(fit), "ggplot")
  fit$preds <- NULL
  expect_error(post_pred_group(fit, grp = "Year"), "add_preds = TRUE")
})

test_that("xlab is passed to the plots", {
  dat <- data_prep(ex_comp, strata, "length", "count")
  fit <- suppressMessages(fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                                    iter = 300, warmup = 100, thin = 2))
  expect_equal(Fx_plot(fit, grp = "Year", xlab = "Age")$labels$x, "Age")
  sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
  expect_equal(scaled_comp_plot(sc, grps = "Year", xlab = "Age")$labels$x, "Age")
  expect_equal(scaled_ridge_plot(sc, grps = "Year", plot_y = "Year", xlab = "Age")$labels$x, "Age")
})

test_that("year effects recover the simulated increase in fish size", {
  dat <- data_prep(ex_comp, strata, "length", "count")
  fit <- suppressMessages(fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                                    iter = 300, warmup = 100, thin = 2))
  sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
  mean_len <- sc |>
    dplyr::group_by(Year, .draw, bin) |>
    dplyr::summarise(t = sum(tot_by_bin), .groups = "drop_last") |>
    dplyr::summarise(m = sum(as.numeric(as.character(bin)) * t) / sum(t), .groups = "drop") |>
    dplyr::group_by(Year) |>
    dplyr::summarise(m = mean(m), .groups = "drop")
  expect_gt(mean_len$m[mean_len$Year == "2020"], mean_len$m[mean_len$Year == "2015"])
})

test_that("spatial workflow runs on the bundled data", {
  dat <- data_prep(ex_comp, c(strata, "x", "y"), "length", "count") |>
    dplyr::mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
  mesh <- suppressMessages(sdmTMB::make_mesh(dat, c("x", "y"), cutoff = 25))
  fit <- suppressMessages(fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                                    mesh = mesh, time = "Year", coords = c("x", "y"),
                                    iter = 300, warmup = 100, thin = 2))
  expect_true(isTRUE(fit$model$sd_report$pdHess))
  sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
  expect_s3_class(sc, "data.frame")
  expect_s3_class(spatial_field_plot(fit, lon_step = 20, lat_step = 20), "ggplot")
})

test_that("scaled_comp_plot handles scaled = FALSE with overlays, and multiple groups", {
  dat <- data_prep(ex_comp, strata, "length", "count")
  fit <- suppressMessages(fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
                                    iter = 300, warmup = 100, thin = 2))
  sc  <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
  raw <- dplyr::rename(ex_comp, bin = length)
  raw$tcatch <- raw$count * 100
  p <- scaled_comp_plot(sc, grps = "Year", scaled = FALSE,
                        comp_are = raw, comp_are2 = raw, cvar = "tcatch")
  expect_s3_class(p, "ggplot")
  # two grouping variables give a facet grid
  p2 <- scaled_comp_plot(sc, grps = c("Year", "Season"))
  expect_s3_class(p2$facet, "FacetGrid")
  # scaled = TRUE with both overlays
  p3 <- scaled_comp_plot(sc, grps = "Year", comp_are = raw, comp_are2 = raw, cvar = "count")
  expect_s3_class(p3, "ggplot")
})
