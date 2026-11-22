# Scenario helpers for the "Scenario Comparison" section of spatial_bias_simulation.Rmd.
#
# A parametrised version of the simulation in that vignette (same truth, selectivity and sampling
# code), used to compare the non-spatial and spatial models under different sampling designs.
# The results used in the vignette are produced by run_scenarios.R and stored in
# scenario_results.csv.
#
# The simulated field is constant over time, so the spatial-only model uses
# fit_model(spatiotemporal = FALSE) (fit_spatial_only()).
#
# Assumptions and settings.
# - Every scenario uses the same seeds (make_world(): 123, make_sample(): 456). Scenarios with
#   different field ranges may therefore have correlated realised fields. Each
#   scenario is a single realisation, so absolute errors are comparable within a scenario and
#   only loosely between scenarios.
# - sdmTMB fits use iter = 400, warmup = 100 and thin = 3, which gives 100 simulation draws in
#   scale_comps(). The draws add Monte Carlo variation to the errors. set.seed(seed_eval) in
#   run_scenario() fixes it for reproducibility.
suppressPackageStartupMessages({library(tidyverse); library(ComPoM); library(sdmTMB)})
options(dplyr.summarise.inform = FALSE)
source("vignette_helpers.R")   # Selectivity()

clo <- function(x) x / sum(x)
n_bins <- 10; M <- 0.2
popfreq <- clo(exp(-(1:n_bins) * M))
L50_base <- 6.5; L95 <- 0.75
SIM_ITER <- 400; SIM_WARMUP <- 100; SIM_THIN <- 3

season_probs_by_year <- function(year, n_years, n_seasons, shift = TRUE) {
  if (!shift) return(rep(1 / n_seasons, n_seasons))
  alpha <- (year - 1) / (n_years - 1)
  w <- numeric(n_seasons)
  w[1] <- 4 * (1 - alpha) + 0.5 * alpha
  w[n_seasons] <- 0.5 * (1 - alpha) + 4 * alpha
  if (n_seasons > 2) w[2:(n_seasons - 1)] <- 0.5
  w / sum(w)
}

# Truth. sigma is the SD of the field on log(L50) and range the range of its exponential
# covariance. trend adds a west-to-east gradient of plus or minus trend across the domain.
# As in the vignette, the catch is drawn independently (gamma, mean 10,000) for every location,
# year, season and bin, and the true catch is the selected abundance times that catch.
make_world <- function(sigma = 0.5, range = 10, n_locations = 250, n_years = 10, n_seasons = 4,
                       seed = 123, trend = 0) {
  set.seed(seed)
  locations <- data.frame(loc_id = 1:n_locations, x = runif(n_locations, 0, 100),
                          y = runif(n_locations, 0, 100))
  d <- as.matrix(dist(cbind(locations$x, locations$y)))
  locations$spatial_effect <- as.numeric(
    MASS::mvrnorm(1, rep(0, n_locations), sigma^2 * exp(-d / range)))
  locations$spatial_effect <- locations$spatial_effect + trend * (locations$x - 50) / 50
  year_effects <- cumsum(rnorm(n_years, 0, 0.1))
  season_effects <- rnorm(n_seasons, 0, 0.2)
  season_effects <- season_effects - mean(season_effects)
  temporal_df <- expand_grid(Year = factor(1:n_years), Season = factor(1:n_seasons)) |>
    mutate(year_effect = year_effects[as.integer(Year)],
           season_effect = season_effects[as.integer(Season)])
  true_catch_by_bin <- expand_grid(loc_id = 1:n_locations, Year = factor(1:n_years),
                                   Season = factor(1:n_seasons), bin = 1:n_bins) |>
    left_join(locations, by = "loc_id") |> left_join(temporal_df, by = c("Year", "Season")) |>
    mutate(L50 = exp(log(L50_base) + spatial_effect + year_effect + season_effect),
           catch = rgamma(n(), 100, 0.01),
           tcatch = Selectivity(bin, L50, L95) * popfreq[bin] * catch)
  list(locations = locations, temporal_df = temporal_df, truth = true_catch_by_bin,
       n_years = n_years, n_seasons = n_seasons)
}

# footprint = "corner" samples locations with probability proportional to exp((x + y) * loc_bias)
# in every year (the vignette design). loc_bias = 0 gives uniform sampling.
# footprint = "moving" uses a Gaussian kernel (SD = width) whose centre moves linearly from c0
# to c1 over the years. seas_shift = FALSE gives equal sampling probability for all seasons.
make_sample <- function(w, footprint = c("corner", "moving"), loc_bias = 1 / 25,
                        c0 = c(85, 85), c1 = c(15, 15), width = 25, seas_shift = TRUE,
                        n_samps = 50, n_per = 100, seed = 456) {
  footprint <- match.arg(footprint)
  set.seed(seed)
  with(w, {
    probs_for <- function(yr) {
      p <- if (footprint == "corner") exp((locations$x + locations$y) * loc_bias) else {
        a <- (yr - 1) / (n_years - 1)
        cc <- (1 - a) * c0 + a * c1
        exp(-((locations$x - cc[1])^2 + (locations$y - cc[2])^2) / (2 * width^2))
      }
      p / sum(p)
    }
    samples <- tibble(samp_id = 1:(n_samps * n_years),
                      Year = factor(rep(1:n_years, each = n_samps))) |>
      mutate(loc_id = unlist(lapply(as.integer(Year), function(yr)
               sample(locations$loc_id, 1, prob = probs_for(yr)))),
             Season = factor(sapply(as.integer(Year), function(yr)
               sample(1:n_seasons, 1,
                      prob = season_probs_by_year(yr, n_years, n_seasons, seas_shift))))) |>
      left_join(locations, by = "loc_id") |> left_join(temporal_df, by = c("Year", "Season")) |>
      mutate(L50 = exp(log(L50_base) + spatial_effect + year_effect + season_effect))
    long <- samples |> group_by(samp_id) |> nest() |>
      mutate(comp_samps = lapply(data, function(df)
        which(rmultinom(n_per, 1, Selectivity(1:n_bins, df$L50, L95) * popfreq) == 1,
              arr.ind = TRUE)[, 1])) |>
      unnest(c(data, comp_samps)) |> mutate(count = 1) |> rename(bin = comp_samps)
    list(long = long)
  })
}

make_mesh_world <- function(w, input_data, cutoff = 3) {
  all_locs_sf <- sf::st_as_sf(w$locations, coords = c("x", "y"))
  bnd <- fmesher::fm_nonconvex_hull(all_locs_sf, convex = -0.1, format = "sf")
  msh <- fmesher::fm_mesh_2d(w$locations[, c("x", "y")], cutoff = cutoff, boundary = bnd)
  make_mesh(input_data, xy_cols = c("x", "y"), mesh = msh)
}

prep <- function(smp, bin_width = NULL) {
  data_prep(smp$long, vars_for_grouping = c("loc_id", "Year", "Season", "x", "y"),
            bin_lab = "bin", bin_width = bin_width) |>
    mutate(x = as.numeric(as.character(x)), y = as.numeric(as.character(y)))
}

# Spatial-only model: fit_model() without the spatio-temporal terms
fit_spatial_only <- function(input_data, mesh, share_sd = TRUE, iter = SIM_ITER,
                             warmup = SIM_WARMUP, thin = SIM_THIN) {
  fit_model(re_form = "Season", ffx_form = "Year", data = input_data, backend = "TMB",
            mesh = mesh, time = "Year", coords = c("x", "y"), spatiotemporal = FALSE,
            share_spatial_sd = share_sd, iter = iter, warmup = warmup, thin = thin)
}

fits_for <- function(w, smp, spatial = c("full", "spatial_only", "none"), bin_width = NULL) {
  input_data <- prep(smp, bin_width)
  out <- list(input_data = input_data)
  out$nonspatial <- fit_model(re_form = "Season", ffx_form = "Year", data = input_data,
                              backend = "TMB", iter = SIM_ITER, warmup = SIM_WARMUP,
                              thin = SIM_THIN)
  if (any(c("full", "spatial_only") %in% spatial)) out$mesh <- make_mesh_world(w, input_data)
  if ("full" %in% spatial)
    out$full <- fit_model(re_form = "Season", ffx_form = "Year", data = input_data,
                          backend = "TMB", mesh = out$mesh, time = "Year",
                          coords = c("x", "y"), iter = SIM_ITER, warmup = SIM_WARMUP,
                          thin = SIM_THIN)
  if ("spatial_only" %in% spatial) out$spatial_only <- fit_spatial_only(input_data, out$mesh)
  out
}

# Error of the annual compositions (relative to the true catch composition) by method.
# MAE is the mean absolute error over years and bins. const is the mean absolute value of the
# per-bin mean error over years. varying is the mean absolute residual after removing that
# per-bin mean. The mean of an absolute value is not additive, so MAE is not const + varying.
# With bin_width, the truth bins are floor(bin / bin_width) * bin_width. This matches the bin
# labels of data_prep(), which also uses the lower edge of the wider bin. data_prep() fills zero
# counts, so every year and bin has a raw proportion.
evaluate <- function(w, smp, fits, label = "", bin_width = NULL) {
  tr <- w$truth
  if (!is.null(bin_width)) tr <- tr |> mutate(bin = floor(bin / bin_width) * bin_width)
  true_year <- tr |> group_by(Year, bin) |> summarise(tcatch = sum(tcatch)) |>
    group_by(Year) |> mutate(truth = tcatch / sum(tcatch)) |> ungroup() |>
    select(Year, bin, truth)
  raw <- fits$input_data |> group_by(Year, bin) |> summarise(v = sum(tot_by_bin)) |>
    group_by(Year) |> mutate(raw = v / sum(v)) |> ungroup() |>
    mutate(bin = as.numeric(as.character(bin))) |> select(Year, bin, raw)
  # The full scale_comps() output has one row per stratum, bin and draw (10 million rows for
  # 10,000 strata, 10 bins and 100 draws). It is therefore computed in 10 chunks of locations
  # and only the annual sums per draw are kept. Every chunk must contain all Years and Seasons,
  # because scale_comps() needs every level of the fixed-effect factor in the scaling data.
  # Each chunk draws its own simulation sample, so draw k of one chunk is independent of draw k
  # of another. Only the mean over draws is used below, which is not affected by this. Interval
  # estimates from these annual sums would be too narrow.
  ann_chunked <- function(fit, sdf, grps, n_chunks = 1) {
    if (n_chunks > 1) {
      key <- paste(sdf$x, sdf$y)
      sdf$.chunk <- as.integer(factor(key)) %% n_chunks
    } else sdf$.chunk <- 0L
    sums <- map_dfr(unique(sdf$.chunk), function(ch) {
      sc <- scale_comps(dplyr::filter(sdf, .chunk == ch) |> select(-.chunk), "catch", fit, grps)
      out <- sc |> group_by(Year, bin, .draw) |> summarise(v = sum(tot_by_bin))
      rm(sc); gc(FALSE)
      out
    }) |> group_by(Year, bin, .draw) |> summarise(v = sum(v))
    sums |> group_by(Year, .draw) |> mutate(p = v / sum(v)) |> group_by(Year, bin) |>
      summarise(p = mean(p)) |> ungroup() |> mutate(bin = as.numeric(as.character(bin)))
  }
  sdf_sp <- tr |> group_by(Season, Year, x, y) |> summarise(catch = sum(tcatch)) |> ungroup()
  # The non-spatial scaling data are restricted to the Year x Season cells with samples, as in
  # the vignette. All 40 cells were sampled in every scenario, so this does not change the
  # population that is scaled. The check below stops a scenario where a cell is unsampled,
  # because the spatial model and the truth cover all cells.
  sdf_ns <- tr |>
    semi_join(fits$input_data |> ungroup() |> distinct(Year, Season), by = c("Year", "Season")) |>
    group_by(Year, Season) |> summarise(catch = sum(tcatch)) |> ungroup()
  stopifnot(nrow(sdf_ns) == dplyr::n_distinct(tr$Year, tr$Season))
  est <- true_year |> left_join(raw, by = c("Year", "bin"))
  est$nonspatial <- (est |> left_join(
    ann_chunked(fits$nonspatial, sdf_ns, c("Year", "Season")), by = c("Year", "bin")))$p
  for (m in intersect(c("full", "spatial_only"), names(fits)))
    est[[m]] <- (est |> left_join(
      ann_chunked(fits[[m]], sdf_sp, c("Season", "Year"), n_chunks = 10),
      by = c("Year", "bin")))$p
  stopifnot(!anyNA(est))   # a missing bin would otherwise propagate NA into the errors
  meths <- setdiff(names(est), c("Year", "bin", "truth"))
  map_dfr(meths, function(m) {
    e <- est[[m]] - est$truth
    ec <- tapply(e, est$bin, mean)
    er <- e - ec[as.character(est$bin)]
    tibble(scenario = label, method = m, MAE = mean(abs(e)), const = mean(abs(ec)),
           varying = mean(abs(er)))
  })
}

# Mean of the simulated field (log L50 scale) over the samples of the first and last year (a
# location sampled twice counts twice), and over all locations of the domain. The Year factor
# levels are in numeric order, because the factor is built from the integers 1 to n_years.
field_at_samples <- function(w, smp) {
  d <- smp$long |> ungroup() |> distinct(samp_id, Year, spatial_effect)
  by_year <- tapply(d$spatial_effect, d$Year, mean)
  c(field_year1 = unname(by_year[1]), field_lastyear = unname(by_year[length(by_year)]),
    field_domain = mean(w$locations$spatial_effect))
}

# -----------------------------------------------------------------------------------------------
# Scenario definitions. The defaults are 250 locations, 10 years, 4 seasons, 50 samples per year
# and 100 fish per sample. footprint is "corner" (north-east in every year) or "moving"
# (north-east to south-west over the years). Uniform sampling has loc_bias = 0 and no seasonal
# shift. range is the range of the simulated field, trend a large-scale west-to-east gradient
# in log(L50), n_per the number of fish per sample and bin_width the aggregation of the bins.
# world and smp may hold further arguments of make_world() and make_sample(). full = TRUE also
# fits the model with spatio-temporal fields, which took about 5 minutes in a test run.
scenario_specs <- list(
  list(name = "A", label = "Vignette design", range = 10, trend = 0, full = TRUE,
       smp = list(footprint = "corner")),
  list(name = "D", label = "Uniform sampling", range = 10, trend = 0,
       smp = list(footprint = "corner", loc_bias = 0, seas_shift = FALSE)),
  list(name = "B", label = "Moving footprint", range = 10, trend = 0,
       smp = list(footprint = "moving")),
  list(name = "E", label = "Longer range (30)", range = 30, trend = 0,
       smp = list(footprint = "corner")),
  list(name = "C", label = "Moving, longer range (30)", range = 30, trend = 0,
       smp = list(footprint = "moving")),
  list(name = "G", label = "300 fish per sample", range = 10, trend = 0,
       smp = list(footprint = "corner", n_per = 300)),
  list(name = "F", label = "Wider bins (2 units)", range = 10, trend = 0,
       smp = list(footprint = "corner"), bin_width = 2),
  list(name = "H1", label = "Gradient", range = 10, trend = 0.4,
       smp = list(footprint = "corner")),
  list(name = "H2", label = "Gradient, moving footprint", range = 10, trend = 0.4,
       smp = list(footprint = "moving")),
  list(name = "H3", label = "Gradient, uniform sampling", range = 10, trend = 0.4,
       smp = list(footprint = "corner", loc_bias = 0, seas_shift = FALSE))
)

# Simulates, fits and evaluates one scenario. Returns one row per method (raw samples,
# non-spatial model and spatial-only model, and for scenario A also the model with
# spatio-temporal fields). sigma_spatial, range_est and the convergence columns describe the
# spatial-only fit and are NA for the other methods.
run_scenario <- function(spec, seed_eval = 1) {
  t0 <- Sys.time()
  w <- do.call(make_world, c(list(range = spec$range, trend = spec$trend), spec$world))
  smp <- do.call(make_sample, c(list(w = w), spec$smp))
  spatial <- if (isTRUE(spec$full)) c("full", "spatial_only") else "spatial_only"
  fits <- suppressWarnings(fits_for(w, smp, spatial = spatial, bin_width = spec$bin_width))
  set.seed(seed_eval)   # simulation draws in scale_comps()
  res <- evaluate(w, smp, fits, label = spec$name, bin_width = spec$bin_width)
  sp_model <- fits$spatial_only$model
  sp_pars <- tidy(sp_model, "ran_pars")
  sp_sanity <- suppressMessages(sanity(sp_model))
  f <- field_at_samples(w, smp)
  is_sp <- res$method == "spatial_only"
  # [1] takes the shared SD (share_sd = TRUE gives one sigma_Z row) and the single range
  res |> mutate(
    label = spec$label,
    n_locs_sampled = dplyr::n_distinct(smp$long$loc_id),
    n_bins = dplyr::n_distinct(fits$input_data$bin),
    sigma_spatial = ifelse(is_sp, sp_pars$estimate[sp_pars$term == "sigma_Z"][1], NA_real_),
    range_est = ifelse(is_sp, sp_pars$estimate[sp_pars$term == "range"][1], NA_real_),
    pd_hessian = ifelse(is_sp, isTRUE(sp_model$sd_report$pdHess), NA),
    sanity_ok = ifelse(is_sp, isTRUE(sp_sanity$all_ok), NA),
    field_year1 = f[["field_year1"]], field_lastyear = f[["field_lastyear"]],
    field_domain = f[["field_domain"]],
    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
