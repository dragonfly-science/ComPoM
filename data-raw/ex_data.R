# Simulated example data for ComPoM (examples, README and tests).
# Run from the package root: source("data-raw/ex_data.R")
set.seed(20261)
lens <- seq(10, 60, by = 10)

strata <- expand.grid(Year = 2015:2020, Season = c("summer", "winter"),
                      Vessel = paste0("V", 1:4), trip = 1:3,
                      stringsAsFactors = FALSE)
strata$x   <- runif(nrow(strata), 0, 100)
strata$y   <- runif(nrow(strata), 0, 100)
strata$sst <- round(runif(nrow(strata), 12, 20), 2)
strata$catch <- round(rlnorm(nrow(strata), log(500), 0.5))

vessel_shift <- c(V1 = -5, V2 = 0, V3 = 5, V4 = 10)   # vessel-specific size selection
comp_for <- function(i) {
  mu <- 30 + 1.5 * (strata$Year[i] - 2015) +             # fish grow larger over the years
        vessel_shift[strata$Vessel[i]] +
        ifelse(strata$Season[i] == "winter", 3, 0) +
        0.3 * (strata$sst[i] - 16) +
        0.12 * (strata$x[i] - 50) - 0.08 * (strata$y[i] - 50)   # spatial size gradient
  p <- dnorm(lens, mu, 12)
  p / sum(p)
}
ex_comp <- do.call(rbind, lapply(seq_len(nrow(strata)), function(i) {
  data.frame(strata[i, c("Year", "Season", "Vessel", "trip", "x", "y", "sst")],
             length = lens, count = as.vector(rmultinom(1, 60, comp_for(i))),
             row.names = NULL)
}))
ex_catch <- strata
ex_comp$Year <- factor(ex_comp$Year); ex_catch$Year <- factor(ex_catch$Year)
for (v in c("Season", "Vessel")) {
  ex_comp[[v]] <- factor(ex_comp[[v]]); ex_catch[[v]] <- factor(ex_catch[[v]])
}
usethis_ok <- requireNamespace("usethis", quietly = TRUE)
save(ex_comp, ex_catch, file = "data/ex_data.rda", compress = "xz", version = 2)
