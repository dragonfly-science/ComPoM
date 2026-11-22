## Regression tests on numerical output. The simulated data have a known true
## composition per year (helper-sim-truth.R); fits and scaled compositions
## must recover it, and the formula strings are fixed.

fit_truth <- suppressMessages(suppressWarnings(
  fit_model("year", re_form = "season", data = sim_truth_dat, backend = "TMB",
            iter = 1000, warmup = 200, thin = 4)
))

# mean absolute difference between the estimated and the true yearly composition
comp_error <- function(sc) {
  est <- sc |>
    dplyr::group_by(year, bin, .draw) |>
    dplyr::summarise(tot = sum(tot_by_bin), .groups = "drop") |>
    dplyr::group_by(year, .draw) |>
    dplyr::mutate(p = tot / sum(tot)) |>
    dplyr::group_by(year, bin) |>
    dplyr::summarise(p = mean(p), .groups = "drop") |>
    dplyr::mutate(bin = as.numeric(as.character(bin)))
  tru <- sim_truth$true_comp
  tru$bin <- as.numeric(tru$bin)
  j <- dplyr::inner_join(est, tru, by = c("year", "bin"), suffix = c("_est", "_true"))
  expect_equal(nrow(j), 4 * 5)
  mean(abs(j$p_est - j$p_true))
}

test_that("TMB catch weighting recovers the true yearly composition", {
  sdf <- expand.grid(year = factor(2001:2004), season = factor(c("S1", "S2", "S3")))
  sdf$catch <- 1000
  sc <- suppressWarnings(scale_comps(sdf, "catch", fit_truth, grps = c("year", "season")))
  expect_lt(comp_error(sc), 0.02)
})

test_that("TMB CPUE scaling recovers the true yearly composition", {
  sdf <- data.frame(year = factor(2001:2004), cpue = 5)
  sc <- suppressWarnings(scale_comps(sdf, "cpue", fit_truth, grps = "year",
                                     remove_nuisance = TRUE))
  expect_lt(comp_error(sc), 0.02)
})

test_that("scaled numbers are proportional to the stratum weight", {
  sdf1 <- data.frame(year = factor(2001:2004), cpue = 1)
  sdf2 <- data.frame(year = factor(2001:2004), cpue = 10)
  set.seed(1)
  a <- suppressWarnings(scale_comps(sdf1, "cpue", fit_truth, grps = "year", remove_nuisance = TRUE))
  set.seed(1)
  b <- suppressWarnings(scale_comps(sdf2, "cpue", fit_truth, grps = "year", remove_nuisance = TRUE))
  expect_equal(b$tot_by_bin, 10 * a$tot_by_bin, tolerance = 1e-8)
})

test_that("predicted numbers sum to approximately the stratum weight", {
  sdf <- data.frame(year = factor(2001:2004), cpue = 100)
  sc <- suppressWarnings(scale_comps(sdf, "cpue", fit_truth, grps = "year",
                                     remove_nuisance = TRUE))
  tot <- sc |> dplyr::group_by(year, .draw) |>
    dplyr::summarise(t = sum(tot_by_bin), .groups = "drop")
  expect_equal(mean(tot$t), 100, tolerance = 0.05)
})

test_that("formula strings are stable", {
  f <- function(...) suppressMessages(parse_form(data = sim_dat, ...))$form
  expect_equal(deparse(f(backend = "TMB", ffx_form = "year", re_form = "season")),
               "tot_by_bin ~ 0 + bin + factor(year):bin + (1 | bin:season)")
  expect_equal(deparse(f(backend = "brms", ffx_form = "year", re_form = "season"))[1],
               "tot_by_bin ~ offset(log(n)) + 0 + bin + factor(year):bin + (1 | ")
  squash <- function(x) gsub("\\s+", "", paste(deparse(x), collapse = ""))
  expect_match(squash(f(backend = "TMB", ffx_form = "year", cvars = "n",
                        spl = "by_factor", knots = "4")),
               "s(n,by=bin,k=4)", fixed = TRUE)
  expect_match(squash(f(backend = "TMB", ffx_form = "year", cvars = "n",
                        spl = "2D", knots = "c(3,3)")),
               "t2(n,cbin,k=c(3,3))", fixed = TRUE)
})

test_that("the legacy 'by' smooth builds one t2 term per non-reference bin", {
  pf <- suppressMessages(parse_form(data = sim_cvar_dat, backend = "TMB", ffx_form = "year",
                                    cvars = "sst", spl = "by", knots = "3"))
  f <- gsub("\\s+", "", paste(deparse(pf$form), collapse = ""))
  expect_equal(lengths(regmatches(f, gregexpr("t2\\(", f))), 3)  # 4 bins -> 3 terms
})
