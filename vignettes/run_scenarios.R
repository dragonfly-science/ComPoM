# Regenerates scenario_results.csv, which the "Scenario Comparison" section of
# spatial_bias_simulation.Rmd reads.
#
# Run from the vignettes directory with Rscript run_scenarios.R. The run starts 2 processes, and
# each fit uses detectCores() %/% 2 TMB threads. The minutes column of scenario_results.csv gives
# the time per scenario. Timings and memory were measured with /usr/bin/time on a machine with
# 20 cores and 62 GB. Scenario A also fits the model with spatio-temporal fields, which took about
# 5 minutes. Peak memory was about 1 GB per process for the other scenarios and 3 GB for
# scenario A, because the scale_comps() output is aggregated in chunks.
#
# SCENARIO_SMOKE=1 Rscript run_scenarios.R runs a small version of scenario A (100 locations,
# 4 years, 3 seasons, 40 samples per year) to check the pipeline. It took about 15 seconds. It
# prints the results and does not write the csv file. It does not cover full = TRUE, bin_width
# or n_per.
source("scenario_helpers.R")

if (identical(Sys.getenv("SCENARIO_SMOKE"), "1")) {
  spec <- scenario_specs[[1]]
  spec$world <- list(n_locations = 100, n_years = 4, n_seasons = 3)
  spec$smp$n_samps <- 40
  spec$full <- FALSE   # the model with spatio-temporal fields takes about 5 minutes
  res <- run_scenario(spec)
  print(as.data.frame(res), digits = 3)
  stopifnot(nrow(res) == 3, all(is.finite(res$MAE)),
            all(res$pd_hessian[res$method == "spatial_only"]))
  quit(save = "no")
}

# mc.preschedule = FALSE starts one process per scenario, two at a time, so that the long
# scenario A does not leave one worker with most of the work
res <- parallel::mclapply(scenario_specs, function(s) try(run_scenario(s)), mc.cores = 2,
                          mc.preschedule = FALSE)
# A worker killed by the operating system returns NULL, and an error returns a try-error
failed <- !vapply(res, is.data.frame, logical(1))
if (any(failed)) {
  names_failed <- vapply(scenario_specs, `[[`, "", "name")[failed]
  print(res[failed])
  stop("scenarios failed: ", paste(names_failed, collapse = ", "))
}
write.csv(dplyr::bind_rows(res), "scenario_results.csv", row.names = FALSE)
