## Tests for post-processing functions
## Fixture objects (sim_fit, sim_scale_df) are created in helper-sim.R

# ── post_pred_group ───────────────────────────────────────────────────────────

test_that("post_pred_group returns a ggplot", {
  p <- suppressWarnings(post_pred_group(sim_fit, grp = "year"))
  expect_s3_class(p, "ggplot")
})

test_that("post_pred_group facets by the supplied group", {
  p <- suppressWarnings(post_pred_group(sim_fit, grp = "year"))
  expect_true(inherits(p$facet, "FacetWrap"))
})

# ── Fx_plot ───────────────────────────────────────────────────────────────────

test_that("Fx_plot returns a ggplot", {
  p <- suppressWarnings(Fx_plot(sim_fit, grp = "year", form = "~(1|bin:year)"))
  expect_s3_class(p, "ggplot")
})

test_that("Fx_plot facets by the supplied group", {
  p <- suppressWarnings(Fx_plot(sim_fit, grp = "year", form = "~(1|bin:year)"))
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("Fx_plot y-axis is labelled 'Multiplier'", {
  p <- suppressWarnings(Fx_plot(sim_fit, grp = "year", form = "~(1|bin:year)"))
  expect_equal(p$labels$y, "Multiplier")
})

# ── scale_comps ───────────────────────────────────────────────────────────────

test_that("scale_comps returns a data frame", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df = sim_scale_df,
      predvar  = "catch",
      fit      = sim_fit,
      grps     = c("year", "area")
    )
  )
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps output contains tot_by_bin and .draw columns", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df = sim_scale_df,
      predvar  = "catch",
      fit      = sim_fit,
      grps     = c("year", "area")
    )
  )
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

test_that("scale_comps produces one row per bin x group x draw", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df = sim_scale_df,
      predvar  = "catch",
      fit      = sim_fit,
      grps     = c("year", "area")
    )
  )
  n_draws  <- sim_fit$nsim
  n_groups <- nrow(sim_scale_df)           # 4 year x area combos
  n_bins   <- length(unique(sim_fit$data$bin))  # 4 bins
  expect_equal(nrow(sc), n_draws * n_groups * n_bins)
})

# ── scaled_comp_plot ──────────────────────────────────────────────────────────

test_that("scaled_comp_plot returns a ggplot (proportions)", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df = sim_scale_df,
      predvar  = "catch",
      fit      = sim_fit,
      grps     = c("year", "area")
    )
  )
  p <- suppressWarnings(
    scaled_comp_plot(scaled_comp = sc, grps = c("year", "area"), scaled = TRUE)
  )
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns a ggplot (counts)", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df = sim_scale_df,
      predvar  = "catch",
      fit      = sim_fit,
      grps     = c("year", "area")
    )
  )
  p <- suppressWarnings(
    scaled_comp_plot(scaled_comp = sc, grps = c("year", "area"), scaled = FALSE)
  )
  expect_s3_class(p, "ggplot")
})
