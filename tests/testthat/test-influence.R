# step_comps(), influ_comps() and their plots, on the two-term fixture (year + area fixed
# effects, season random effect) from helper-sim.R
inf_grps <- c("year", "area", "season")
inf_grid <- data.frame(year = factor(2020:2021))

test_that("comp_quantile and the stat argument", {
  expect_equal(comp_quantile(0.5)(c(10, 20, 30), c(0.2, 0.5, 0.3)), 20)
  expect_equal(comp_quantile(0.9)(c(10, 20, 30), c(0.2, 0.5, 0.3)), 30)
  expect_error(step_comps(sim_fit_two, sim_scale_df, "catch", inf_grps, "year", stat = "nope"),
               "stat")
})

test_that("step_comps adds terms in formula order and the first step is flat", {
  st <- step_comps(sim_fit_two, scale_df = sim_scale_df, predvar = "catch", grps = inf_grps,
                   by = "year", grid = inf_grid, selectivity = c("area", "season"))
  expect_s3_class(st, "comp_step")
  expect_equal(levels(st$steps$step), c("bin only", "+ year", "+ area", "+ season"))
  expect_setequal(unique(st$steps$component), c("catch", "index"))
  flat <- dplyr::filter(st$steps, step == "bin only", component == "index")
  expect_equal(diff(flat$mean), 0, tolerance = 1e-8)
  expect_equal(nrow(st$summary), 4)
  expect_true(all(c("df", "logLik", "AIC") %in% names(st$summary)))
})

test_that("step_comps honours `order`, `stat` and a single component", {
  st <- step_comps(sim_fit_two, scale_df = sim_scale_df, predvar = "catch", grps = inf_grps,
                   by = "year", order = c("season", "area", "year"), stat = 0.9)
  expect_equal(levels(st$steps$step), c("bin only", "+ season", "+ area", "+ year"))
  expect_equal(unique(st$steps$component), "catch")
  expect_match(st$stat_label, "90th")
  expect_error(step_comps(sim_fit_two, scale_df = sim_scale_df, grps = inf_grps, by = "year",
                          order = c("year", "area")), "order")
})

test_that("step_comps refuses the brms backend", {
  skip_if_no_brms()
  expect_error(step_comps(sim_fit_brms, scale_df = sim_scale_df, grps = inf_grps, by = "year"),
               "TMB")
})

test_that("influ_comps returns bin-resolved influence and zero index influence for selectivity terms", {
  inf <- influ_comps(sim_fit_two, scale_df = sim_scale_df, predvar = "catch", grps = inf_grps,
                     by = "year", grid = inf_grid, selectivity = c("area", "season"))
  expect_s3_class(inf, "comp_influ")
  expect_setequal(names(inf$groups), c("area", "season"))
  expect_true(all(is.finite(inf$influence$mean)))
  idx <- dplyr::filter(inf$influence, component == "index")
  expect_true(all(idx$mean == 0))
  cat_inf <- dplyr::filter(inf$influence, component == "catch", term == "area")
  expect_gt(max(abs(cat_inf$mean)), 0)
  expect_true(all(c("var", "level", "bin", "effect", "lo", "hi") %in% names(inf$effects)))
  expect_true(all(c("term", "share") %in% names(inf$distribution)))
  expect_error(influ_comps(sim_fit_two, scale_df = sim_scale_df, grps = inf_grps, by = "year",
                           terms = "year"), "cannot be a term")
})

test_that("plot_step_comps and plot_cdi_comps build", {
  st <- step_comps(sim_fit_two, scale_df = sim_scale_df, predvar = "catch", grps = inf_grps,
                   by = "year", grid = inf_grid, selectivity = c("area", "season"))
  expect_no_error(ggplot2::ggplot_build(plot_step_comps(st)))
  inf <- influ_comps(sim_fit_two, scale_df = sim_scale_df, predvar = "catch", grps = inf_grps,
                     by = "year", grid = inf_grid, selectivity = c("area", "season"))
  p <- plot_cdi_comps(inf, "area")
  expect_s3_class(p, "patchwork")
  pdf(NULL); on.exit(dev.off())
  expect_no_error(print(p))
  expect_error(plot_cdi_comps(inf, "nope"), "must be one of")
})
