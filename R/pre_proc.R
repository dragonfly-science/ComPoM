#' Prepare compositional count data
#'
#' Takes a dataset of counts (number of sampled individuals) per bin
#' (length/age/etc) across a range of strata (\code{vars_for_grouping}) and
#' expands it for use in the Poisson-multinomial model. Counts are summed
#' within stratum and bin, bins without observations within a stratum are
#' filled with zero counts, and the stratum total \code{n} (used as offset in
#' \code{fit_model()}) is added.
#'
#' Strata are defined by the unique combinations of \code{vars_for_grouping}.
#' Any variable later used in \code{fit_model()} (fixed effects, nuisance
#' random effects, continuous covariates, coordinates, time) must be included
#' here, otherwise it is aggregated away.
#' @param comp A compositional dataset with columns for counts (can be names differently) per bin (can be named differently), and all `vars_for_grouping`
#' @param vars_for_grouping A character vector of variables used for grouping the compositional dataframe
#' @param bin_lab A string with the column label for the bins if not "bin"
#' @param count_var A string with the column label for the counts by bins if not "count"
#' @param min_n Integer. Minimum total count per grouping stratum. Strata with
#'   fewer than \code{min_n} total fish are dropped before fitting. Useful for
#'   removing sparse early years that would create all-zero design matrix
#'   columns and a singular Hessian. Default \code{1L} (only empty strata dropped).
#' @param bin_width Numeric. If supplied, bins are aggregated into wider bins of
#'   this width before summarising (e.g. \code{bin_width = 20} merges 10cm bins
#'   into 20cm bins). Applied before \code{min_n} filtering. Default \code{NULL}
#'   (no aggregation).
#' @return A grouped data frame with one row per stratum and bin, containing
#'   \code{vars_for_grouping}, \code{bin}, the count per bin
#'   (\code{tot_by_bin}) and the stratum total (\code{n}).
#' @export
#' @import dplyr
#' @import tidyr
#' @import rlang

data_prep <- function(
    comp = NULL,
    vars_for_grouping = NULL,
    bin_lab          = "bin",
    count_var        = 'count',
    min_n            = 1L,
    bin_width        = NULL
){

  counts <- sym(count_var)
  bin    <- sym(bin_lab)

  # Optionally aggregate to wider bins
  if (!is.null(bin_width)) {
    comp <- comp |>
      dplyr::mutate(!!bin := floor(!!bin / bin_width) * bin_width)
  }

  ## bin counts by grouping vars (0 counts are included via complete())
  comp |>
    group_by(across(all_of(c(vars_for_grouping, bin_lab)))) |>
    summarize(tot_by_bin = sum(!!counts, na.rm = TRUE), .groups = "drop") |>
    complete(nesting(!!!syms(vars_for_grouping)), !!bin,
             fill = list(tot_by_bin = 0)) |>
    group_by(across(all_of(vars_for_grouping))) |>
    mutate(n = sum(tot_by_bin)) |>
    rename(bin = bin_lab) |>
    filter(n >= min_n)
}

#' Parse a formula
#'
#' Constructs the Poisson-multinomial model formula, separating fixed effects
#' (variables of primary interest, entered as \code{factor(var):bin}
#' interactions) from nuisance random effects (entered as \code{(1|bin:var)}
#' IID random intercepts). Called internally by \code{fit_model()}; it is
#' exported to allow inspection of the formula and model data before fitting.
#'
#' The linear predictor always contains bin-specific intercepts
#' (\code{0 + bin}). Continuous covariates in \code{cvars} are added as smooths
#' according to \code{spl}:
#' \describe{
#'   \item{\code{"2D"}}{\code{t2(cvar, cbin, k = knots)}, a tensor-product
#'     smooth over the covariate and the numeric bin value. Requires numeric
#'     bins (e.g. length).}
#'   \item{\code{"by_factor"}}{\code{s(cvar, by = bin, k = knots)}, a separate
#'     1D smooth of the covariate for each bin. Also applicable to
#'     non-ordered categories (e.g. species).}
#'   \item{\code{"by"}}{Legacy option. Separate \code{t2()} smooths for each
#'     non-reference bin via bin dummy columns.}
#' }
#'
#' @param data Prepared data frame from \code{data_prep()}
#' @param backend the backend used for model fitting (\code{"brms"} or \code{"TMB"})
#' @param ffx_form Character string of primary fixed-effect variables separated
#'   by \code{+}, e.g. \code{"Year"}. These become \code{factor(Year):bin}
#'   fixed-effect interactions in the model.
#' @param re_form Character string of nuisance random-effect variables separated
#'   by \code{+}, e.g. \code{"Season + Gear"}. These become
#'   \code{(1|bin:Season)} IID random intercepts. \code{scale_comps()}
#'   evaluates them at each stratum's level (\code{remove_nuisance = FALSE})
#'   or removes them (\code{remove_nuisance = TRUE}).
#' @param cvars Continuous variable to be fitted as 2D spline with bin
#' @param knots knots (1D or 2D) for cvar spline; as character for now
#' @param spl 1D spline over cvar per bin, or 2D spline over bins and cvar; as character for now
#' @param coords Character vector of coordinate column names to preserve (TMB backend only)
#' @param time Time column name to preserve (TMB backend only)
#' @return A list with the model formula (\code{form}), the model data with
#'   grouping variables converted to factors (\code{data}), and the names of
#'   the nuisance random-effect variables (\code{re_vars}).
#' @import stringr
#' @import recipes
#' @export
#'
parse_form <- function(data,
                       backend  = 'brms',
                       ffx_form = 'Year',
                       re_form  = NULL,
                       cvars    = NULL,
                       spl      = "2D",
                       knots    = 'c(5,10)',
                       coords   = NULL,
                       time     = NULL
                       ){

  # ── fixed-effect part: factor(var):bin interactions ─────────────────────────
  ffx_terms <- NULL
  if (!is.null(ffx_form) && nchar(trimws(ffx_form)) > 0) {
    ffx_parts <- stringr::str_remove_all(
      stringr::str_split(ffx_form, '\\+')[[1]], pattern = ' '
    )
    ffx_terms <- paste0('factor(', ffx_parts, '):bin')
  }

  # ── nuisance random-effect part: (1|bin:var) IID intercepts ─────────────────
  re_terms <- NULL
  re_vars  <- NULL
  if (!is.null(re_form) && nchar(trimws(re_form)) > 0) {
    re_vars  <- stringr::str_remove_all(
      stringr::str_split(re_form, '\\+')[[1]], pattern = ' '
    )
    re_terms <- paste0('(1|bin:', re_vars, ')')
  }

  newform <- paste(
    c('0 + bin', ffx_terms, re_terms),
    collapse = ' + '
  )

  if (is.numeric(data$bin)) data$bin <- factor(data$bin, levels = sort(unique(data$bin)))

  if (!is.null(cvars) & spl == '2D') {
    data$cbin <- as.numeric(as.character(data$bin))
    for (i in 1:length(cvars))
      newform <- paste(newform, paste("t2(", cvars[i], ", cbin, k = ", knots, ")"), sep = ' + ')
  }
  if (!is.null(cvars) & spl == 'by') {
    mm      <- model.matrix(data = data, ~bin)
    binvars <- unique(data$bin)
    colnames(mm) <- binvars
    lb <- length(binvars)
    data <- cbind(data, mm[, 2:lb])
    for (i in 1:length(cvars))
      for (b in binvars[2:lb])
        newform <- paste(newform,
                         paste("t2(", cvars[i], sprintf(", by=%s, k = ", b), knots, ")"),
                         sep = ' + ')
  }
  if (!is.null(cvars) & spl == 'by_factor') {
    # s(cvar, by=bin, k=knots): separate 1D smooth of cvar for each bin level.
    # Full-rank alongside factor(Year):bin interactions since sst varies
    # within year across SST strata, independently of the year dummies.
    for (i in seq_along(cvars))
      newform <- paste(newform,
                       paste0("s(", cvars[i], ", by=bin, k=", knots, ")"),
                       sep = ' + ')
  }

  form <- switch(backend,
    brms = as.formula(paste("tot_by_bin~offset(log(n)) +", newform)),
    TMB  = as.formula(paste("tot_by_bin~", newform))
  )

  print(form)

  # ── factor conversion: preserve coords, time, cvars ─────────────────────────
  exclude_pattern <- paste(c('.*_by_.*', '^n$', cvars, if (!is.null(cvars)) '^cbin$'),
                           collapse = '|')
  if (!is.null(coords))
    exclude_pattern <- paste(exclude_pattern,
                             paste0('^', coords, '$', collapse = '|'), sep = '|')
  if (!is.null(time))
    exclude_pattern <- paste(exclude_pattern, paste0('^', time, '$'), sep = '|')

  data <- data |> dplyr::ungroup() |>
    dplyr::mutate(dplyr::across(-dplyr::matches(exclude_pattern), as.factor))

  if (backend == 'TMB') {
    # Build interaction columns for IID RE terms
    if (!is.null(re_vars)) {
      for (v in re_vars) {
        col <- paste0('bin:', v)
        data[[col]] <- paste(data$bin, data[[v]], sep = '_')
      }
    }

    # Also build factor(var):bin columns (already handled by formula dummy coding,
    # but we need the factor cols present)
    exclude_pattern_tmb <- paste(c('.*_by_.*', '^n$', cvars, if (!is.null(cvars)) '^cbin$'),
                                 collapse = '|')
    if (!is.null(coords))
      exclude_pattern_tmb <- paste(exclude_pattern_tmb,
                                   paste0('^', coords, '$', collapse = '|'), sep = '|')
    if (!is.null(time))
      exclude_pattern_tmb <- paste(exclude_pattern_tmb, paste0('^', time, '$'), sep = '|')

    data <- data |>
      dplyr::mutate(dplyr::across(-dplyr::matches(exclude_pattern_tmb), as.factor))
  }

  list(
    form    = form,
    data    = data,
    re_vars = re_vars   # nuisance RE variable names, used by scale_comps
  )
}

#' Fit the Poisson-multinomial model
#'
#' Builds the model formula with \code{parse_form()} and fits it with either
#' \code{brms} (Bayesian, \code{backend = "brms"}) or \code{sdmTMB} (maximum
#' likelihood, \code{backend = "TMB"}). Counts per stratum and bin are modelled
#' as Poisson with \code{log(n)} as offset, so that the implied within-stratum
#' proportions follow a multinomial-logit model (see
#' \code{vignette("ComPoM")} for the model equations).
#'
#' With \code{backend = "TMB"} and a \code{mesh}, the model adds one spatial
#' Gaussian random field per bin and one spatio-temporal field per bin and
#' time step (set up with \code{sdmTMB::make_category_svc()}).
#'
#' With \code{backend = "brms"}, fixed effects receive a
#' \code{student_t(10, 0, 1)} prior. With \code{backend = "TMB"}, uncertainty
#' in downstream functions is obtained from \code{nsim} draws from the
#' approximate joint sampling distribution of the parameters.
#'
#' @param ffx_form A character string of primary fixed-effect variables (RHS),
#'   e.g. \code{"Year"}. Passed to \code{parse_form()}.
#' @param re_form A character string of nuisance random-effect variables,
#'   e.g. \code{"Season + Gear"}. These become \code{(1|bin:var)} IID
#'   intercepts. \code{NULL} = none.
#' @param data Prepared data frame from \code{data_prep()}.
#' @param dist Response distribution passed to \code{brms::brm()} as
#'   \code{family} (brms backend only; the TMB backend always uses
#'   \code{poisson()}).
#' @param backend the backend used for model fitting: \code{"brms"} or
#'   \code{"TMB"}.
#' @param cvars Continuous covariate(s) to be fitted as smooths with bin; see
#'   \code{spl}.
#' @param knots Basis dimension(s) for the \code{cvars} smooth, as a character
#'   string (e.g. \code{"4"} or \code{"c(5,10)"} for \code{spl = "2D"}).
#' @param spl Smooth type for \code{cvars}: \code{"2D"}, \code{"by_factor"} or
#'   \code{"by"}; see \code{parse_form()}.
#' @param chains,cores,threads Number of MCMC chains, cores and within-chain
#'   threads (brms backend only).
#' @param brms_backend Stan interface used by brms (\code{"rstan"} or
#'   \code{"cmdstanr"}).
#' @param iter,warmup,thin MCMC iterations, warm-up iterations and thinning
#'   (brms). For the TMB backend these only set the number of simulation draws
#'   used downstream, \code{nsim = (iter - warmup) / thin}.
#' @param mtd,ad Stan \code{max_treedepth} and \code{adapt_delta} (brms only).
#' @param refresh Stan progress refresh interval (brms only).
#' @param add_preds Add posterior predictive draws for the fitted data,
#'   required by \code{post_pred_group()}. Can be slow.
#' @param mesh An sdmTMB mesh object created with make_mesh() (TMB backend only)
#' @param time Column name for time (year) variable (TMB backend only)
#' @param share_spatial_sd Share spatial SD across bins (TMB backend only)
#' @param share_spatiotemporal_sd Share spatiotemporal SD across bins (TMB backend only)
#' @param coords Character vector of coordinate column names (TMB backend only)
#' @param priors An \code{sdmTMBpriors()} object for regularisation (TMB
#'   non-spatial backend only; currently not passed to the spatial fit).
#'   Useful when the design matrix is sparse (e.g. \code{normal(0, 2)} on fixed effects).
#' @param newton_loops Number of extra Newton optimisation steps after initial fit
#'   (TMB non-spatial only). Default 1. Increase to 2-3 for difficult convergence.
#'   For spatial fits use \code{sdmTMB::run_extra_optimization()} on
#'   \code{$model}.
#' @return An object of class \code{ComPoM}: a list with the fitted model
#'   (\code{model}, a \code{brmsfit} or \code{sdmTMB} object), the model data
#'   (\code{data}), the number of draws used downstream (\code{nsim}),
#'   posterior predictive draws (\code{preds}, or \code{NULL}), and the
#'   settings needed by the post-processing functions (\code{mesh},
#'   \code{spatial_varying}, \code{time}, \code{coords}, \code{re_vars},
#'   \code{cvars}).
#' @seealso \code{data_prep()}, \code{scale_comps()}, \code{Fx_plot()},
#'   \code{vignette("ComPoM")}
#' @export
#' @import brms
#' @import sdmTMB
#' @import tidybayes
#'
fit_model <- function(ffx_form,
                      re_form  = NULL,
                      cvars    = NULL,
                      spl      = '2D',
                      knots    = 'c(5,10)',
                      data     = NULL,
                      dist     = 'poisson()',
                      chains   = 4,
                      threads  = 16,
                      cores    = 32,
                      backend  = 'brms',
                      brms_backend = 'rstan',
                      iter     = 1000,
                      warmup   = 500,
                      mtd      = 10,
                      ad       = 0.8,
                      thin     = 2,
                      refresh  = 10,
                      add_preds = FALSE,
                      mesh     = NULL,
                      time     = NULL,
                      share_spatial_sd = TRUE,
                      share_spatiotemporal_sd = TRUE,
                      coords   = c("x", "y"),
                      priors   = NULL,
                      newton_loops = 1) {

  spatial_varying <- !is.null(mesh)

  df <- parse_form(
    data     = data,
    backend  = backend,
    ffx_form = ffx_form,
    re_form  = re_form,
    cvars    = cvars,
    spl      = spl,
    knots    = knots,
    coords   = if (backend == 'TMB' && spatial_varying) coords else NULL,
    time     = if (backend == 'TMB') time else NULL
  )

  form    <- df$form
  data    <- df$data
  re_vars <- df$re_vars   # nuisance RE variable names (may be NULL)

  fixed <- paste0("bin", unique(data$bin)[1])
  pp    <- set_prior("student_t(10, 0, 1)", class = "b")

  if (backend == 'brms') {
    mod <- brm(bf(form),
               family       = dist,
               prior        = pp,
               data         = data,
               backend      = brms_backend,
               chains       = chains,
               threads      = threads,
               cores        = cores,
               iter         = iter,
               warmup       = warmup,
               thin         = thin,
               refresh      = refresh,
               control      = list(max_treedepth = mtd, adapt_delta = ad)
    )
    cat('Converged: ', all(brms::rhat(mod) < 1.05))

  } else if (backend == 'TMB') {

    if (spatial_varying) {
      if (is.null(time)) stop("'time' must be provided for spatial models")

      svc_setup <- sdmTMB::make_category_svc(
        data                    = data,
        category_column         = "bin",
        time_column             = time,
        share_spatial_sd        = share_spatial_sd,
        share_spatiotemporal_sd = share_spatiotemporal_sd
      )

      control_sdmTMB <- sdmTMBcontrol(
        map          = svc_setup$svc_map,
        parallel     = parallel::detectCores(),
        multiphase   = TRUE,
        newton_loops = 0,
        #profile      = "b_j",
        getsd        = TRUE
      )

      mod <- sdmTMB(
        data             = svc_setup$data_expanded,
        offset           = log(data$n),
        formula          = form,
        family           = poisson(),
        mesh             = mesh,
        spatial          = "off",
        spatiotemporal   = "off",
        time             = time,
        spatial_varying  = svc_setup$svc_formula,
        control          = control_sdmTMB
      )
    } else {
      mod <- sdmTMB(
        data    = data,
        offset  = log(data$n),
        formula = form,
        family  = poisson(),
        spatial      = "off",
        priors       = if (!is.null(priors)) priors else sdmTMBpriors(),
        control      = sdmTMBcontrol(
          newton_loops = newton_loops#,
          #profile      = "b_j"
        )
      )
    }
  }

  if (add_preds)
    pp_draws <- switch(backend,
      brms = add_predicted_draws(data, mod),
      TMB  = simulate(mod,
                      nsim             = (iter - warmup) / thin,
                      type             = "mle-mvn",
                      mle_mvn_samples  = "multiple"
               ) |>
               as.data.frame() |>
               cbind(data) |>
               dplyr::mutate(.row = dplyr::row_number()) |>
               tidyr::pivot_longer(
                 cols      = dplyr::matches("V[0-9]+"),
                 names_to  = "iter",
                 values_to = ".prediction"
               ) |>
               dplyr::mutate(.draw = as.numeric(gsub("V", "", iter)))
    )
  else pp_draws <- NULL

  out <- list(
    model           = mod,
    data            = data,
    nsim            = (iter - warmup) / thin,
    preds           = pp_draws,
    mesh            = mesh,
    spatial_varying = spatial_varying,
    time            = time,
    coords          = if (spatial_varying) coords else NULL,
    re_vars         = re_vars,   # stored for use in scale_comps
    cvars           = cvars      # continuous spline covariates
  )

  class(out) <- 'ComPoM'
  attr(out, 'call') <- match.call()
  out
}



