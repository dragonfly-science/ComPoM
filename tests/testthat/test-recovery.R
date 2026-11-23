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
                                     mode = "index", weight = "cpue"))
  expect_lt(comp_error(sc), 0.02)
})

test_that("scaled numbers are proportional to the stratum weight", {
  sdf1 <- data.frame(year = factor(2001:2004), cpue = 1)
  sdf2 <- data.frame(year = factor(2001:2004), cpue = 10)
  set.seed(1)
  a <- suppressWarnings(scale_comps(sdf1, "cpue", fit_truth, grps = "year", mode = "index", weight = "cpue"))
  set.seed(1)
  b <- suppressWarnings(scale_comps(sdf2, "cpue", fit_truth, grps = "year", mode = "index", weight = "cpue"))
  expect_equal(b$tot_by_bin, 10 * a$tot_by_bin, tolerance = 1e-8)
})

test_that("predicted numbers sum to approximately the stratum weight", {
  sdf <- data.frame(year = factor(2001:2004), cpue = 100)
  sc <- suppressWarnings(scale_comps(sdf, "cpue", fit_truth, grps = "year",
                                     mode = "index", weight = "cpue"))
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

test_that("index mode with a selectivity fixed effect recovers the reference-gear composition", {
  set.seed(11)
  bins <- c(10, 20, 30, 40, 50)
  # gear B selects larger fish than gear A; the population composition shifts with year
  comp_for <- function(yr, gear) {
    p <- dnorm(bins, mean = 25 + 3 * yr + if (gear == "B") 10 else 0, sd = 9)
    p / sum(p)
  }
  strata <- expand.grid(year = factor(1:3), gear = c("A", "B"), rep = 1:6)
  samp <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i) {
    p <- comp_for(as.integer(strata$year[i]), as.character(strata$gear[i]))
    data.frame(strata[rep(i, length(bins)), ], bin = bins,
               count = as.vector(rmultinom(1, 400, p)), unit = i, row.names = NULL)
  }))
  dat <- data_prep(samp, c("year", "gear", "unit"), "bin", "count")
  fit <- suppressMessages(suppressWarnings(
    fit_model("year + gear", data = dat, backend = "TMB", iter = 300, warmup = 100, thin = 1)
  ))
  sc <- suppressWarnings(scale_comps(fit = fit, mode = "index", grps = "year",
                                     grid = data.frame(year = factor(1:3)),
                                     selectivity = "gear", reference = list(gear = "A")))
  est <- sc |>
    dplyr::group_by(year, bin, .draw) |>
    dplyr::summarise(v = sum(tot_by_bin), .groups = "drop") |>
    dplyr::group_by(year, .draw) |>
    dplyr::mutate(p = v / sum(v)) |>
    dplyr::group_by(year, bin) |>
    dplyr::summarise(p = mean(p), .groups = "drop")
  truth <- do.call(rbind, lapply(1:3, function(y)
    data.frame(year = factor(y), bin = factor(bins), pt = comp_for(y, "A"))))
  est$bin <- factor(as.numeric(as.character(est$bin)))
  cmp <- dplyr::inner_join(est, truth, by = c("year", "bin"))
  expect_equal(nrow(cmp), 15)
  expect_lt(mean(abs(cmp$p - cmp$pt)), 0.02)
})

test_that("spatial index mode recovers the population composition under spatially biased sampling", {
  set.seed(5)
  bins <- c(10, 20, 30, 40, 50)
  # mean length increases with x and year; gear B selects larger fish
  comp_for <- function(yr, x, gear) {
    p <- dnorm(bins, 25 + 3 * yr + 0.15 * (x - 50) + if (gear == "B") 10 else 0, 9)
    p / sum(p)
  }
  # sampled locations cluster in the east; the population grid is regular
  locs   <- data.frame(loc = 1:80, x = 100 * rbeta(80, 1.5, 1), y = runif(80, 0, 100))
  strata <- merge(expand.grid(loc = 1:80, year = factor(1:3), gear = c("A", "B")), locs, by = "loc")
  samp <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i) {
    p <- comp_for(as.integer(strata$year[i]), strata$x[i], as.character(strata$gear[i]))
    data.frame(strata[rep(i, length(bins)), ], bin = bins,
               count = as.vector(rmultinom(1, 300, p)), unit = i, row.names = NULL)
  }))
  dat <- data_prep(samp, c("unit", "year", "gear", "x", "y"), "bin", "count") |>
    dplyr::mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
  mesh <- suppressMessages(sdmTMB::make_mesh(dat, c("x", "y"), cutoff = 20))
  fit <- suppressMessages(suppressWarnings(
    fit_model("year + gear", data = dat, backend = "TMB", mesh = mesh, time = "year",
              coords = c("x", "y"), spatiotemporal = FALSE, iter = 300, warmup = 100, thin = 1)
  ))
  fit$model <- suppressMessages(suppressWarnings(
    sdmTMB::run_extra_optimization(fit$model, newton_loops = 2)
  ))

  cells <- expand.grid(x = seq(5, 95, 10), y = seq(5, 95, 10), year = factor(1:3))
  sc <- suppressWarnings(scale_comps(fit = fit, mode = "index", grps = "year", grid = cells,
                                     selectivity = "gear", reference = list(gear = "A")))
  est <- sc |>
    dplyr::group_by(year, bin, .draw) |>
    dplyr::summarise(v = sum(tot_by_bin), .groups = "drop") |>
    dplyr::group_by(year, .draw) |>
    dplyr::mutate(p = v / sum(v)) |>
    dplyr::group_by(year, bin) |>
    dplyr::summarise(p = mean(p), .groups = "drop") |>
    dplyr::mutate(bin = as.numeric(as.character(bin)))
  # true population composition: mean over the regular grid, reference gear A
  truth <- do.call(rbind, lapply(1:3, function(y) {
    data.frame(year = factor(y), bin = bins,
               pt = rowMeans(sapply(seq_len(100), function(i) comp_for(y, cells$x[i], "A"))))
  }))
  raw <- samp[samp$gear == "A", ] |>
    dplyr::group_by(year, bin) |>
    dplyr::summarise(v = sum(count), .groups = "drop") |>
    dplyr::group_by(year) |>
    dplyr::mutate(p = v / sum(v)) |>
    dplyr::ungroup()

  mae_index <- mean(abs(dplyr::inner_join(est, truth, by = c("year", "bin")) |>
                          dplyr::mutate(e = p - pt) |> dplyr::pull(e)))
  mae_raw   <- mean(abs(dplyr::inner_join(raw, truth, by = c("year", "bin")) |>
                          dplyr::mutate(e = p - pt) |> dplyr::pull(e)))
  expect_lt(mae_index, 0.015)
  # the raw samples overrepresent the east, so the index composition is closer to the truth
  expect_lt(mae_index, 0.6 * mae_raw)
})
