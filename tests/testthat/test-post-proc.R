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
      grps     = c("year", "area", "season")
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# post_pred_group
# ══════════════════════════════════════════════════════════════════════════════

test_that("post_pred_group returns a ggplot [brms]", {
  skip_if_no_brms()
  p <- suppressWarnings(post_pred_group(sim_fit_brms, grp = "year"))
  expect_s3_class(p, "ggplot")
})

test_that("post_pred_group facets by group [brms]", {
  skip_if_no_brms()
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
  skip_if_no_brms()
  p <- suppressWarnings(Fx_plot(sim_fit_brms, grp = "year", form = "~(1|bin:year)"))
  expect_s3_class(p, "ggplot")
})

test_that("Fx_plot y-axis is 'Multiplier' [brms]", {
  skip_if_no_brms()
  p <- suppressWarnings(Fx_plot(sim_fit_brms, grp = "year", form = "~(1|bin:year)"))
  expect_equal(p$labels$y, "Multiplier")
})

test_that("Fx_plot works for a non-fixed-effect group [brms]", {
  skip_if_no_brms()
  # year is held at its mode, leaving factor(year) with one level in newdata
  p <- suppressWarnings(Fx_plot(sim_fit_brms, grp = "season"))
  expect_s3_class(p, "ggplot")
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
  skip_if_no_brms()
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps contains expected columns [brms]", {
  skip_if_no_brms()
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

# ── scale_comps index mode ───────────────────────────────────────────────

# Predicted proportion per draw x stratum x bin, for comparing seasons
season_props <- function(sc) {
  sc |>
    dplyr::ungroup() |>
    dplyr::mutate(p = tot_by_bin / n) |>
    dplyr::select(year, area, season, bin, .draw, p) |>
    tidyr::pivot_wider(names_from = season, values_from = p)
}

test_that("catch weighting keeps each stratum's nuisance level [TMB non-spatial]", {
  sc <- get_sc(sim_fit_tmb, sim_scale_df)
  expect_setequal(as.character(unique(sc$season)), c("S1", "S2"))
})

test_that("catch weighting errors when a nuisance variable is missing [TMB]", {
  expect_error(
    scale_comps(sim_scale_df, predvar = "catch", fit = sim_fit_tmb,
                grps = c("year", "area")),
    "season"
  )
})

test_that("catch weighting errors when a nuisance variable is missing [brms]", {
  skip_if_no_brms()
  expect_error(
    scale_comps(sim_scale_df, predvar = "catch", fit = sim_fit_brms,
                grps = c("year", "area")),
    "season"
  )
})

test_that("catch weighting warns on unseen nuisance levels [TMB non-spatial]", {
  sdf <- sim_scale_df
  levels(sdf$season) <- c("S1", "S9")
  expect_warning(
    scale_comps(sdf, predvar = "catch", fit = sim_fit_tmb,
                grps = c("year", "area", "season")),
    "S9"
  )
})

test_that("index mode gives identical compositions across nuisance levels [TMB non-spatial]", {
  sc <- suppressWarnings(
    scale_comps(sim_scale_df, predvar = "catch", fit = sim_fit_tmb,
                grps = c("year", "area", "season"), mode = "index", weight = "catch")
  )
  sp <- season_props(sc)
  expect_equal(sp$S1, sp$S2)
})

test_that("index mode keeps the stratum columns in the output [TMB non-spatial]", {
  sc <- suppressWarnings(
    scale_comps(sim_scale_df, predvar = "catch", fit = sim_fit_tmb,
                grps = c("year", "area", "season"), mode = "index", weight = "catch")
  )
  expect_setequal(as.character(unique(sc$season)), c("S1", "S2"))
})

test_that("index mode runs [brms]", {
  skip_if_no_brms()
  sc <- suppressWarnings(
    scale_comps(sim_scale_df, predvar = "catch", fit = sim_fit_brms,
                grps = c("year", "area"), mode = "index", weight = "catch")
  )
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

test_that("scale_comps index mode returns a data frame [TMB non-spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df,
      predvar         = "catch",
      fit             = sim_fit_tmb,
      grps            = c("year", "area"),
      mode = "index", weight = "catch"
    )
  )
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps index mode contains expected columns [TMB non-spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df,
      predvar         = "catch",
      fit             = sim_fit_tmb,
      grps            = c("year", "area"),
      mode = "index", weight = "catch"
    )
  )
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

test_that("scale_comps index mode returns a data frame [TMB spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df_spatial,
      predvar         = "catch",
      fit             = sim_fit_tmb_spatial,
      grps            = c("year", "area"),
      mode = "index", weight = "catch"
    )
  )
  expect_s3_class(sc, "data.frame")
})

test_that("scale_comps index mode contains expected columns [TMB spatial]", {
  sc <- suppressWarnings(
    scale_comps(
      scale_df        = sim_scale_df_spatial,
      predvar         = "catch",
      fit             = sim_fit_tmb_spatial,
      grps            = c("year", "area"),
      mode = "index", weight = "catch"
    )
  )
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})

# ══════════════════════════════════════════════════════════════════════════════

test_that("scaled_comp_plot returns ggplot — proportions [brms]", {
  skip_if_no_brms()
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(scaled_comp_plot(sc, grps = c("year", "area"), scaled = TRUE))
  expect_s3_class(p, "ggplot")
})

test_that("scaled_comp_plot returns ggplot — counts [brms]", {
  skip_if_no_brms()
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
  skip_if_no_brms()
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year")
  )
  expect_s3_class(p, "ggplot")
})

test_that("scaled_ridge_plot errors without plot_y [brms]", {
  skip_if_no_brms()
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  expect_error(
    scaled_ridge_plot(sc, grps = c("year", "area")),
    "`plot_y` must be supplied"
  )
})

test_that("scaled_ridge_plot facets when plot_facet is supplied [brms]", {
  skip_if_no_brms()
  sc <- get_sc(sim_fit_brms, sim_scale_df)
  p  <- suppressWarnings(
    scaled_ridge_plot(sc, grps = c("year", "area"), plot_y = "year", plot_facet = "area")
  )
  expect_s3_class(p, "ggplot")
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("scaled_ridge_plot bin_range subsets bins [brms]", {
  skip_if_no_brms()
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

test_that("theme_compom and legend labels in scaled_comp_plot", {
  expect_s3_class(theme_compom(), "theme")
  sc <- scale_comps(sim_scale_df, "catch", sim_fit_tmb, grps = c("year", "season"))
  ref <- dplyr::rename(sim_comp, tcatch = count)
  p <- scaled_comp_plot(sc, grps = "year", comp_are = ref, comp_label = "True")
  b <- ggplot2::ggplot_build(p)
  expect_true("True" %in% b$plot$scales$get_scales("colour")$get_labels())
})

test_that("post_pred_group plot builds", {
  expect_no_error(ggplot2::ggplot_build(post_pred_group(sim_fit_tmb_spatial, grp = "year")))
})

# ── index mode: roles, references, deprecated alias ───────────────────────────

test_that("remove_nuisance = TRUE still works and matches mode = 'index'", {
  set.seed(1)
  a <- suppressWarnings(scale_comps(sim_scale_df, "catch", sim_fit_tmb,
                                    grps = c("year", "area"), remove_nuisance = TRUE))
  set.seed(1)
  b <- suppressWarnings(scale_comps(sim_scale_df, "catch", sim_fit_tmb, grps = c("year", "area"),
                                    mode = "index", weight = "catch"))
  expect_equal(a$tot_by_bin, b$tot_by_bin)
})

test_that("index mode uses equal weights when `weight` is NULL", {
  g <- data.frame(year = factor(2020:2021), area = factor("A"))
  sc <- suppressWarnings(scale_comps(fit = sim_fit_tmb, mode = "index", grps = c("year", "area"),
                                     grid = g))
  expect_true(all(sc$n == 1))
})

test_that("index mode errors on terms without a role", {
  g <- data.frame(year = factor(2020:2021))
  expect_error(
    scale_comps(fit = sim_fit_two, mode = "index", grps = "year", grid = g),
    "no role"
  )
})

test_that("index mode checks selectivity and density arguments", {
  g <- data.frame(year = factor(2020:2021), area = factor("A"))
  expect_error(scale_comps(fit = sim_fit_tmb, mode = "index", grps = c("year", "area"),
                           grid = g, selectivity = "year"), "must not be in `grps`")
  expect_error(scale_comps(fit = sim_fit_tmb, mode = "index", grps = c("year", "area"),
                           grid = g, selectivity = "nope"), "not in the model")
  expect_error(scale_comps(fit = sim_fit_tmb, mode = "index", grps = c("year", "area"),
                           grid = g, density = "nope"), "columns of `grid`")
  expect_error(scale_comps(fit = sim_fit_tmb, mode = "index", grps = c("year", "area"),
                           grid = g, weight = "nope"), "not in `grid`")
})

test_that("a selectivity fixed effect is set to the reference and ignored in the grid", {
  g <- data.frame(year = factor(2020:2021))
  set.seed(1)
  a <- suppressWarnings(scale_comps(fit = sim_fit_two, mode = "index", grps = "year", grid = g,
                                    selectivity = "area", reference = list(area = "A")))
  set.seed(1)
  b <- suppressWarnings(scale_comps(fit = sim_fit_two, mode = "index", grps = "year", grid = g,
                                    selectivity = "area", reference = list(area = "B")))
  expect_s3_class(a, "data.frame")
  # different references give different predictions when area has an effect
  expect_false(isTRUE(all.equal(a$tot_by_bin, b$tot_by_bin)))
})

test_that("keeping all random effects via `density` evaluates them at grid values", {
  g <- data.frame(year = factor(2020:2021), area = factor("A"),
                  season = factor(c("S1", "S2")))
  sc <- suppressWarnings(scale_comps(fit = sim_fit_tmb, mode = "index",
                                     grps = c("year", "area", "season"), grid = g,
                                     density = "season"))
  expect_setequal(as.character(unique(sc$season)), c("S1", "S2"))
})

test_that("index mode runs for a spatial grid [TMB spatial]", {
  sc <- suppressWarnings(scale_comps(sim_scale_df_spatial, fit = sim_fit_tmb_spatial,
                                     mode = "index", grps = c("year", "area"),
                                     weight = "catch"))
  expect_true(all(c("tot_by_bin", ".draw", "bin", "x", "y") %in% names(sc)))
})

test_that("index mode runs with a selectivity random effect and a reference [brms]", {
  skip_if_no_brms()
  g <- data.frame(year = factor(2020:2021), area = factor("A"))
  sc <- suppressWarnings(scale_comps(fit = sim_fit_brms, mode = "index",
                                     grps = c("year", "area"), grid = g))
  expect_true(all(c("tot_by_bin", ".draw", "bin") %in% names(sc)))
})
