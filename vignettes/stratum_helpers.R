# Helpers for the stratum-intercept comparison vignette (stratum_intercept.Rmd).
# Not part of the ComPoM package.
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2); library(sdmTMB)})

# Simulate multinomial length samples from strata (year x vessel x trip). The composition shifts
# towards larger bins over the years, and vessels shift the composition up or down. `dispersion`
# multiplies each stratum's proportions by gamma noise (bin-level overdispersion).
sim_strata <- function(seed = 1, nyear = 6, nves = 12, ntrip = 6, bins = seq(20, 90, by = 10),
                       sd_ves = 0.5, nmean = 40, dispersion = 0) {
  set.seed(seed)
  mu_year <- 45 + seq(0, 15, length.out = nyear)
  ves_eff <- rnorm(nves, 0, sd_ves)
  st <- expand.grid(Year = seq_len(nyear), Vessel = seq_len(nves), trip = seq_len(ntrip))
  st <- st[runif(nrow(st)) < 0.6, ]
  st$stratum <- seq_len(nrow(st))
  st$n <- 1 + rpois(nrow(st), nmean)
  p_true <- function(m, v = 0) {
    w <- dnorm(bins, m + 8 * v, 12)
    w / sum(w)
  }
  rows <- lapply(seq_len(nrow(st)), function(i) {
    p <- p_true(mu_year[st$Year[i]], ves_eff[st$Vessel[i]])
    if (dispersion > 0) {
      p <- p * rgamma(length(p), 1 / dispersion, 1 / dispersion)
      p <- p / sum(p)
    }
    data.frame(st[i, c("Year", "Vessel", "trip", "stratum")], bin = bins,
               count = as.vector(rmultinom(1, st$n[i], p)))
  })
  # truth: composition at the average vessel, by year (columns), bins in rows
  list(comp = bind_rows(rows), truth = sapply(mu_year, p_true), bins = bins)
}

# Prepare the data as factors, with one 0/1 column per year (not the first) x bin combination.
# `ref_bin = TRUE` leaves out the first bin, which is needed when a free stratum intercept is added.
make_data <- function(dat, ref_bin = FALSE) {
  d <- as.data.frame(dat)
  d <- d |> mutate(across(-c(tot_by_bin, n, bin), as.factor), bin = as.factor(bin))
  if ("Vessel" %in% names(d)) d$`bin:Vessel` <- interaction(d$bin, d$Vessel, sep = "_")
  bins <- levels(d$bin)
  yrs  <- levels(d$Year)
  cols <- character()
  for (y in yrs[-1]) for (b in (if (ref_bin) bins[-1] else bins)) {
    nm <- paste0("y", y, "_b", b)
    d[[nm]] <- as.numeric(d$Year == y & d$bin == b)
    cols <- c(cols, nm)
  }
  attr(d, "ycols") <- cols
  d
}

# Fit one formulation with sdmTMB. variant: "offset" (the ComPoM model), "random" (random stratum
# intercept), "fixed_naive" (free stratum intercept, all year x bin effects) and "fixed_ref" (free
# stratum intercept, year x bin effects for all bins but the first).
fit_variant <- function(dat, variant = c("offset", "random", "fixed_naive", "fixed_ref"),
                        nuisance = TRUE, family = poisson()) {
  variant <- match.arg(variant)
  d <- make_data(dat, ref_bin = variant == "fixed_ref")
  rhs <- paste("0 + bin +", paste0("`", attr(d, "ycols"), "`", collapse = " + "))
  if (nuisance) rhs <- paste(rhs, "+ (1 | bin:Vessel)")
  rhs <- switch(variant,
    offset = rhs, random = paste(rhs, "+ (1 | stratum)"),
    fixed_naive = , fixed_ref = paste(rhs, "+ stratum"))
  form <- stats::as.formula(paste("tot_by_bin ~", rhs))
  err <- NULL
  tm <- system.time(
    m <- tryCatch(suppressWarnings(sdmTMB(form, data = d, offset = log(d$n), family = family,
                                          spatial = "off")),
                  error = function(e) {err <<- conditionMessage(e); NULL}))[["elapsed"]]
  list(model = m, data = d, variant = variant, nuisance = nuisance, time = tm, error = err,
       pd = !is.null(m) && isTRUE(m$sd_report$pdHess),
       npar = if (is.null(m)) NA else length(m$tmb_obj$par),
       loglik = if (is.null(m)) NA else as.numeric(stats::logLik(m)))
}

# Year composition at the average vessel (random effects off) and a fixed reference stratum.
year_props <- function(fit) {
  if (is.null(fit$model)) return(NULL)
  d <- fit$data
  bins <- levels(d$bin); yrs <- levels(d$Year)
  nd <- expand.grid(bin = factor(bins, levels = bins), Year = factor(yrs, levels = yrs))
  for (cn in attr(d, "ycols")) nd[[cn]] <- 0
  for (y in yrs[-1]) for (b in bins) {
    cn <- paste0("y", y, "_b", b)
    if (cn %in% attr(d, "ycols")) nd[[cn]] <- as.numeric(nd$Year == y & nd$bin == b)
  }
  nd$stratum <- factor(levels(d$stratum)[1], levels = levels(d$stratum))
  nd$Vessel  <- factor(levels(d$Vessel)[1], levels = levels(d$Vessel))
  nd$`bin:Vessel` <- factor(paste(nd$bin, nd$Vessel, sep = "_"), levels = levels(d$`bin:Vessel`))
  eta <- predict(fit$model, newdata = nd, re_form = NA, re_form_iid = NA,
                 offset = rep(0, nrow(nd)))$est
  nd$p <- exp(eta)
  nd <- nd |> group_by(Year) |> mutate(p = p / sum(p)) |> ungroup()
  matrix(nd$p, nrow = length(bins))
}

# Standard deviation of the stratum random intercept (NA if there is none)
stratum_sd <- function(fit) {
  if (is.null(fit$model) || fit$variant != "random") return(NA_real_)
  tt <- tidy(fit$model, "ran_pars")
  tt$estimate[tt$group_name == "stratum"]
}

# Dispersion parameter of a negative binomial fit (NA otherwise)
nb_phi <- function(fit) {
  if (is.null(fit$model) || !grepl("nbinom", fit$model$family$family)) return(NA_real_)
  exp(fit$model$sd_report$par.fixed[["ln_phi"]])
}

# Annual composition from the fixed effects alone (fields and random effects at zero)
coef_props <- function(model, d) {
  b <- tidy(model, "fixed")
  est <- setNames(b$estimate, b$term)
  bins <- levels(d$bin); yrs <- levels(d$Year)
  eta <- sapply(yrs, function(y) sapply(bins, function(bb) {
    e <- est[[paste0("bin", bb)]]
    cn <- paste0("y", y, "_b", bb)
    if (cn %in% names(est)) e <- e + est[[cn]]
    e
  }))
  apply(exp(eta), 2, function(x) x / sum(x))
}

# Simulate multinomial samples at random locations. A Gaussian random field (exponential
# covariance, `range`, SD `sd_field`) shifts the mean of the composition at each location.
sim_spatial <- function(seed = 1, nyear = 4, n_loc = 40, bins = seq(20, 70, by = 10),
                        range = 25, sd_field = 0.6, nmean = 50) {
  set.seed(seed)
  loc <- data.frame(x = runif(n_loc, 0, 100), y = runif(n_loc, 0, 100))
  D <- as.matrix(dist(loc))
  loc$u <- as.numeric(MASS::mvrnorm(1, rep(0, n_loc), sd_field^2 * exp(-D / range)))
  mu_year <- 40 + 4 * (seq_len(nyear) - 1)
  p_of <- function(m, u = 0) { w <- dnorm(bins, m + 8 * u, 12); w / sum(w) }
  # a stratum is a sampling event at a location in a year (new events every year)
  st <- expand.grid(Year = seq_len(nyear), loc = seq_len(n_loc), rep = 1:2)
  st <- st[runif(nrow(st)) < 0.7, ]
  st$stratum <- seq_len(nrow(st)); st$x <- loc$x[st$loc]; st$y <- loc$y[st$loc]
  st$n <- 1 + rpois(nrow(st), nmean)
  rows <- lapply(seq_len(nrow(st)), function(i)
    data.frame(st[i, c("Year", "stratum", "x", "y")], bin = bins,
               count = as.vector(rmultinom(1, st$n[i], p_of(mu_year[st$Year[i]], loc$u[st$loc[i]])))))
  list(comp = bind_rows(rows), truth = sapply(mu_year, p_of), bins = bins)
}

# Spatial-only fit (one spatial field per bin, as fit_model(spatiotemporal = FALSE)) of the
# offset, random and free-intercept formulations.
fit_spatial_variant <- function(dat, mesh, variant = c("offset", "random", "fixed_ref")) {
  variant <- match.arg(variant)
  d <- make_data(dat, ref_bin = variant == "fixed_ref")
  d$x <- as.numeric(as.character(d$x)); d$y <- as.numeric(as.character(d$y))
  svc <- sdmTMB::make_category_svc(data = d, category_column = "bin", time_column = "Year",
                                   share_spatial_sd = TRUE, share_spatiotemporal_sd = TRUE)
  sp_cols <- svc$info$spatial_terms
  svc_formula <- stats::as.formula(paste0("~ 0 + `", paste(sp_cols, collapse = "` + `"), "`"))
  svc_map <- list(ln_tau_Z = factor(rep(1, length(sp_cols))))
  rhs <- paste("0 + bin +", paste0("`", attr(d, "ycols"), "`", collapse = " + "))
  rhs <- switch(variant, offset = rhs, random = paste(rhs, "+ (1 | stratum)"),
                fixed_ref = paste(rhs, "+ stratum"))
  err <- NULL
  tm <- system.time(m <- tryCatch(suppressWarnings(sdmTMB(
    stats::as.formula(paste("tot_by_bin ~", rhs)), data = svc$data_expanded, offset = log(d$n),
    family = poisson(), mesh = mesh, spatial = "off", spatiotemporal = "off", time = "Year",
    spatial_varying = svc_formula,
    control = sdmTMBcontrol(map = svc_map, multiphase = TRUE, newton_loops = 0, getsd = TRUE))),
    error = function(e) {err <<- conditionMessage(e); NULL}))[["elapsed"]]
  list(model = m, data = d, variant = variant, time = tm, error = err,
       pd = !is.null(m) && isTRUE(m$sd_report$pdHess),
       npar = if (is.null(m)) NA else length(m$tmb_obj$par),
       loglik = if (is.null(m)) NA else as.numeric(stats::logLik(m)))
}
