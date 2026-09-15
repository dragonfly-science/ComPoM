## Tests for post-processing functions
## Fixture objects are created in helper-sim.R
##   sim_fit_brms        — brms backend
##   sim_fit_tmb         — sdmTMB non-spatial backend
##   sim_fit_tmb_spatial — sdmTMB spatial backend

# ── Helpers ───────────────────────────────────────────────────────────────────

# Compute scale_comps once per backend to avoid repetition in each test
get_sc <- function(fit, scale_df) {
  suppressWarnings(
    scale_comps(
      scale_df = scale_df,
      predvar  = "catch",
      fit      = fit,
      grps     = c("year", "area")
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# post_pred_group
# ══════════════════════════════════════════════════════════════════════════════

test_that("post_pred_group returns a ggplot [brms]", {
  p <- suppressWarnings(post_pred_group(sim_fit_brms, grp = "year"))
  expect_s3_class(p, "ggplot")
})

test_that("post_pred_group facets by group [brms]", {
  p <- suppressWarnings(post_pred_group(sim_fit_brms, grp = "year"))
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("post_pred_group returns a ggplot [TMB non-spatial]", {
  p <- suppressWarnings(post_pred_group(sim_fit_tmb, grp = "year"))
  expect_s3_class(p, "ggplot")
})

test_that("post_pred_group facets by group [TMB non-spatial]", {
  p <- suppressWarnings(post_pred_group(sim_fit_tmb, grp = "year"))
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("post_pred_group returns a ggplot [TMB spatial]", {
  p <- suppressWarnings(post_pred_group(sim_fit_tmb_spatial, grp = "year"))
  expect_s3_class(p, "ggplot")
})

test_that("post_pred_group facets by group [TMB spatial]", {
  p <- suppressWarnings(post_pred_group(sim_fit_tmb_spatial, grp = "year"))
  expect_true(inherits(p$facet, "FacetWrap"))
})

# ══════════════════════════════════════════════════════════════════════════════
# Fx_plot
# ══════════════════════════════════════════════════════════════════════════════

test_that("Fx_plot returns a ggplot [brms]", {
  p <- suppressWarnings(Fx_plot(sim_fit_brms, grp = "year", form = "~(1|bin:year)"))
  expect_s3_class(p, "ggplot")
})

test_that("Fx_plot y-axis is 'Multiplier' [brms]", {
  p <- suppressWarnings(Fx_plot(sim_fit_brms, grp = "year", form = "~(1|bin:year)"))
  expect_equal(p$labels$y, "Multiplier")
})

test_that("Fx_plot returns a ggplot [TMB non-spatial]", {
  p <- suppressWarnings(Fx_plot(sim_fit_tmb, grp = "year", form = "~(1|bin:year)"))
  expect_s3_class(p, "ggplot")
})

test_that("Fx_plot y-axis is 'Multiplier' [TMB non-spatial]", {
  p <- suppressWarnings(Fx_plot(sim_fit_tmb, grp = "year", form = "~(1|bin:year)"))
  expect_equal(p$labels$y, "Multiplier")
})

test_that("Fx_plot returns a ggplot [TMB spatial]", {
  p <- suppressWarnings(Fx_plot(sim_fit_tmb_spatial, grp = "year", form = "~(1|bin:year)"))
  expect_s3_class(p, "ggplot")
})

test_that("Fx_plot y-axis is 'Multiplier' [TMB spatial]", {
  p <- suppressWarnings(Fx_plot(sim_fit_tmb_spatial, grp = "year", form = "~(1|bin:year)"))
  expect_equal(p$labels$y, "Multiplier")
})

# ══════════════════════════════════════════════════════════════════════════════
# scale_comps
# ══════════════════════════════════════════════════════════════════════════════

test_that("scale_comps returns a data frame [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps contains expected columns [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

test_that("scale_comps returns a data frame [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps contains expected columns [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

test_that("scale_comps correct number of rows [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  n_draws  <- sim_fit_tmb$nsim
  n_groups <- nrow(sim_scale_df)
  n_bins   <- length(unique(sim_fit_tmb$data$bin))
  expect_equal(nrow(sc), n_draws * n_groups * n_bins)
})

test_that("scale_comps returns a data frame [TMB spatial]", {
  sc <- get_sc(sim_fit_tmb_spatial, sim_scale_df_spatial)
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps contains expected columns [TMB spatial]", {
  sc <- get_sc(sim_fit_tmb_spatial, sim_scale_df_spatial)
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

# ── scale_comps remove_nuisance ───────────────────────────────────────────────

test_that("scale_comps remove_nuisance returns a data frame [TMB non-spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df,
      predvar         = "catch",
      fit             = sim_fit_tmb,
      grps            = c("year", "area"),
      remove_nuisance = TRUE
    )
  )
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps remove_nuisance contains expected columns [TMB non-spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df,
      predvar         = "catch",
      fit             = sim_fit_tmb,
      grps            = c("year", "area"),
      remove_nuisance = TRUE
    )
  )
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

test_that("scale_comps remove_nuisance returns a data frame [TMB spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df_spatial,
      predvar         = "catch",
      fit             = sim_fit_tmb_spatial,
      grps            = c("year", "area"),
      remove_nuisance = TRUE
    )
  )
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps remove_nuisance contains expected columns [TMB spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df_spatial,
      predvar         = "catch",
      fit             = sim_fit_tmb_spatial,
      grps            = c("year", "area"),
      remove_nuisance = TRUE
    )
  )
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

# ══════════════════════════════════════════════════════════════════════════════

test_that("scaled_comp_plot returns ggplot — proportions [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns ggplot — counts [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns ggplot — proportions [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns ggplot — counts [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = FALSE))
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns ggplot — proportions [TMB spatial]", {
  sc <- get_sc(sim_fit_tmb_spatial, sim_scale_df_spatial)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns ggplot — counts [TMB spatial]", {
  sc <- get_sc(sim_fit_tmb_spatial, sim_scale_df_spatial)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = FALSE))
  expect_s3_class(p, "ggplot")
})

# ══════════════════════════════════════════════════════════════════════════════
# scaled_ridge_plot
# ══════════════════════════════════════════════════════════════════════════════

test_that("scaled_ridge_plot returns a ggplot [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year")
  )
  expect_s3_class(p, "ggplot")
})

test_that("scaled_ridge_plot errors without plot_y [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  expect_error(
    scaled_ridge_plot(sc, grps = c("year", "area")),
    "`plot_y` must be supplied"
  )
})

test_that("scaled_ridge_plot facets when plot_facet is supplied [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year", plot_facet = "area")
  )
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("scaled_ridge_plot bin_range subsets bins [brms]", {
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year", bin_range = c(20, 40))
  )
  built <- ggplot2::ggplot_build(p)
  bins_plotted <- unique(built$data[[1]]$x)
  expect_true(all(bins_plotted >= 20 & bins_plotted <= 40))
})

test_that("scaled_ridge_plot returns a ggplot [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year")
  )
  expect_s3_class(p, "ggplot")
})

test_that("scaled_ridge_plot returns a ggplot [TMB spatial]", {
  sc <- get_sc(sim_fit_tmb_spatial, sim_scale_df_spatial)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year")
  )
  expect_s3_class(p, "ggplot")
})
