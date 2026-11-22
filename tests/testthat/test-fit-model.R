## Tests for fit_model(): structure of the returned object and input guards.

test_that("fit_model returns a ComPoM object with the documented components", {
  for (f in Filter(Negate(is.null), list(sim_fit_brms, sim_fit_tmb, sim_fit_tmb_spatial))) {
    expect_s3_class(f, "ComPoM")
    expect_true(all(c("model", "data", "nsim", "preds", "mesh", "spatial_varying",
                      "time", "coords", "re_vars", "cvars") %in% names(f)))
  }
})

test_that("fit_model backends return the expected model classes", {
  skip_if_no_brms()
  expect_s3_class(sim_fit_brms$model, "brmsfit")
  expect_s3_class(sim_fit_tmb$model, "sdmTMB")
  expect_s3_class(sim_fit_tmb_spatial$model, "sdmTMB")
})

test_that("spatial flag, coords and nsim are stored", {
  expect_false(sim_fit_tmb$spatial_varying)
  expect_true(sim_fit_tmb_spatial$spatial_varying)
  expect_equal(sim_fit_tmb_spatial$coords, c("x", "y"))
  expect_equal(sim_fit_tmb$nsim, 200)
  expect_null(sim_fit_tmb$coords)
})

test_that("add_preds = FALSE gives no predictive draws", {
  f <- suppressMessages(suppressWarnings(
    fit_model("year", data = sim_dat, backend = "TMB", iter = 400, warmup = 200, thin = 1)
  ))
  expect_null(f$preds)
})

test_that("spatial models require a time variable", {
  expect_error(
    suppressMessages(fit_model("year", data = sim_dat_spatial, backend = "TMB",
                               mesh = sim_mesh, coords = c("x", "y"))),
    "'time' must be provided"
  )
})

test_that("cvars are stored on the fit", {
  expect_equal(sim_fit_cvar$cvars, "sst")
})

test_that("default cores do not exceed chains or available cores", {
  cores <- eval(formals(fit_model)$cores, list(chains = 4))
  expect_lte(cores, 4)
  expect_gte(cores, 1)
})

test_that("the default brms family is accepted by brms", {
  dist <- formals(fit_model)$dist
  expect_no_error(brms::brmsfamily(dist))
  # the family object form documented for `dist` must also work
  expect_no_error(brms::brmsfamily(poisson()$family))
})

test_that("fit_model warns when the sdmTMB Hessian is not positive definite", {
  # 4 locations x 2 years with one spatio-temporal field per bin is not identifiable
  # (this is the former spatial test fixture, which gave a non-positive-definite Hessian)
  set.seed(42)
  d <- expand.grid(year = factor(2020:2021), area = factor(c("A", "B")),
                   season = factor(c("S1", "S2")), bin = c(10, 20, 30, 40))
  d$count <- rpois(nrow(d), lambda = 10 + as.integer(d$bin) / 5)
  d$x <- ifelse(d$area == "A", 0, 10)
  d$y <- ifelse(d$year == "2020", 0, 10)
  dat <- data_prep(d, c("year", "area", "season", "x", "y"), "bin", "count") |>
    dplyr::mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
  mesh <- suppressMessages(sdmTMB::make_mesh(dat, c("x", "y"), cutoff = 3))
  # sdmTMB raises its own convergence warnings; let only ComPoM's reach expect_warning()
  keep_ours <- function(w)
    if (!grepl("The Hessian of the sdmTMB fit", conditionMessage(w))) invokeRestart("muffleWarning")
  expect_warning(
    withCallingHandlers(
      suppressMessages(fit_model("year", re_form = "season", data = dat, backend = "TMB",
                                 mesh = mesh, time = "year", coords = c("x", "y"),
                                 iter = 20, warmup = 10, thin = 1)),
      warning = keep_ours),
    "not positive definite"
  )
})

test_that("[TMB] spatiotemporal = FALSE fits spatial fields only", {
  fit_s <- suppressMessages(suppressWarnings(
    fit_model(ffx_form = "year", re_form = "season", data = sim_dat_spatial,
              backend = "TMB", iter = 400, warmup = 200, thin = 1,
              mesh = sim_mesh, time = "year", coords = c("x", "y"),
              spatiotemporal = FALSE)
  ))
  n_bins <- length(unique(sim_dat_spatial$bin))
  n_time <- length(unique(sim_dat_spatial$year))
  expect_s3_class(fit_s, "ComPoM")
  expect_true(fit_s$spatial_varying)
  # z_i holds the spatially varying terms: one per bin, plus one per bin and year in the full model
  expect_equal(ncol(fit_s$model$tmb_data$z_i), n_bins)
  expect_equal(ncol(sim_fit_tmb_spatial$model$tmb_data$z_i), n_bins * (1 + n_time))
  sc <- suppressWarnings(scale_comps(sim_scale_df_spatial, "catch", fit_s,
                                     grps = c("year", "season", "x", "y")))
  expect_true(all(is.finite(sc$tot_by_bin)))
})
