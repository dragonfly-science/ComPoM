## Covariate smooths and spatial field plots

test_that("cvar_smooth_plot returns a ggplot", {
  p <- suppressWarnings(cvar_smooth_plot(sim_fit_cvar, "sst"))
  expect_s3_class(p, "ggplot")
  expect_gt(nrow(p$data), 0)
})

test_that("cvar_smooth_plot defaults to the first covariate", {
  p <- suppressWarnings(cvar_smooth_plot(sim_fit_cvar))
  expect_s3_class(p, "ggplot")
})

test_that("cvar_smooth_plot errors for a fit without covariates", {
  expect_error(cvar_smooth_plot(sim_fit_tmb), "No continuous covariate")
})

test_that("scale_comps for a covariate model reflects the covariate", {
  # the covariate is averaged within each stratum (it is not a grouping variable),
  # so contrasts in sst are expressed through separate strata
  sdf <- rbind(
    data.frame(year = factor(2001:2003), scen = "low",  sst = 11, catch = 100),
    data.frame(year = factor(2001:2003), scen = "high", sst = 19, catch = 100)
  )
  sc <- suppressWarnings(scale_comps(sdf, "catch", sim_fit_cvar, grps = c("year", "scen")))
  expect_s3_class(sc, "data.frame")
  prop <- sc |>
    dplyr::group_by(year, scen, bin, .draw) |>
    dplyr::summarise(tot = sum(tot_by_bin), .groups = "drop") |>
    dplyr::group_by(year, scen, .draw) |>
    dplyr::mutate(p = tot / sum(tot)) |>
    dplyr::group_by(scen, bin) |>
    dplyr::summarise(p = mean(p), .groups = "drop")
  # the simulated truth has the largest bin more common at high sst
  p40 <- function(s) prop$p[prop$scen == s & as.character(prop$bin) == "40"]
  expect_gt(p40("high"), p40("low"))
})

test_that("spatial_field_plot returns a ggplot for a spatial fit", {
  p <- suppressWarnings(spatial_field_plot(sim_fit_tmb_spatial, lon_step = 20, lat_step = 20))
  expect_s3_class(p, "ggplot")
})

test_that("spatial_field_plot errors for a non-spatial fit", {
  expect_error(spatial_field_plot(sim_fit_tmb), "spatial sdmTMB")
})

# ── 2D (tensor-product) smooth over covariate and bin ─────────────────────────
# In simulations (6 bins, 60-450 strata, knots c(4,4)) the sdmTMB fit of the 2D smooth had a
# non-positive-definite Hessian, so simulated draws (and plots/scaling that need them) fail.
# Only the model set-up is tested here. See NEWS / issue tracker.

test_that("a 2D smooth over covariate and bin can be set up and fitted", {
  d <- expand.grid(year = factor(2001:2003), unit = 1:30, bin = seq(10, 60, by = 10))
  set.seed(7)
  d$sst   <- rep(runif(90, 10, 20), times = 6)
  d$count <- rpois(nrow(d), exp(3 + 0.15 * (d$sst - 15) * (as.integer(factor(d$bin)) - 3.5)))
  dat <- data_prep(d, c("year", "unit", "sst"), "bin", "count") |>
    dplyr::mutate(sst = as.numeric(as.character(sst)))
  pf <- suppressMessages(parse_form(dat, backend = "TMB", ffx_form = "year",
                                    cvars = "sst", spl = "2D", knots = "c(4,4)"))
  expect_true("cbin" %in% names(pf$data))
  fit <- suppressMessages(suppressWarnings(
    fit_model("year", cvars = "sst", spl = "2D", knots = "c(4,4)", data = dat,
              backend = "TMB", iter = 40, warmup = 20, thin = 1)))
  expect_s3_class(fit, "ComPoM")
  expect_equal(fit$cvars, "sst")
})
