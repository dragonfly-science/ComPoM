## Fixtures with a known true composition, for recovery / regression tests.
## Counts are large so that the maximum-likelihood fit is close to the truth.

sim_truth <- local({
  set.seed(2026)
  bins  <- c(10, 20, 30, 40, 50)
  years <- factor(2001:2004)
  # true composition per year: location of a discretised normal shifts with year
  true_comp <- do.call(rbind, lapply(seq_along(years), function(i) {
    p <- dnorm(bins, mean = 25 + 4 * (i - 1), sd = 10)
    data.frame(year = years[i], bin = bins, p = p / sum(p))
  }))
  # 6 sampled strata per year (season x replicate), 300 fish each
  strata <- expand.grid(year = years, season = factor(c("S1", "S2", "S3")), rep = 1:2)
  samp <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i) {
    tc <- true_comp[true_comp$year == strata$year[i], ]
    cbind(strata[i, c("year", "season", "rep")], bin = bins,
          count = as.vector(rmultinom(1, 300, tc$p)))
  }))
  samp$unit <- factor(interaction(samp$season, samp$rep))
  list(bins = bins, true_comp = true_comp, samp = samp)
})

sim_truth_dat <- data_prep(sim_truth$samp, c("year", "season", "unit"), "bin", "count")

# Covariate fixture: one continuous covariate (sst) that shifts the composition
sim_cvar <- local({
  set.seed(7)
  d <- expand.grid(year = factor(2001:2003), unit = 1:20, bin = c(10, 20, 30, 40))
  d$sst   <- rep(runif(60, 10, 20), times = 4)
  d$b     <- as.integer(factor(d$bin))
  d$count <- rpois(nrow(d), exp(3 + 0.15 * (d$sst - 15) * (d$b - 2.5)))
  d
})
sim_cvar_dat <- data_prep(sim_cvar, c("year", "unit", "sst"), "bin", "count") |>
  dplyr::mutate(sst = as.numeric(as.character(sst)))

sim_fit_cvar <- suppressMessages(suppressWarnings(
  fit_model(ffx_form = "year", cvars = "sst", spl = "by_factor", knots = "4",
            data = sim_cvar_dat, backend = "TMB",
            iter = 400, warmup = 200, thin = 1)
))
