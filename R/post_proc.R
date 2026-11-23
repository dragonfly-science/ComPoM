#' Posterior predictive check by group
#'
#' Compares observed counts per bin with posterior predictive (brms) or
#' simulated (sdmTMB) counts, summed within levels of \code{grp}. Predicted
#' totals are rescaled to the observed total within each group and draw, so
#' the check is on the composition rather than on the total.
#' Requires a fit from \code{fit_model(..., add_preds = TRUE)}.
#'
#' @param mod An object from fit_model
#' @param grp Name(s) of the grouping variable(s) in the model data to check
#'   over (e.g. \code{"Year"}). \code{NULL} pools all strata.
#' @param xlab x-axis label (the bin variable)
#' @return A ggplot with predicted medians and intervals (black) and observed
#'   counts (blue), faceted by \code{grp}.
#'
#' @seealso \code{fit_model()}
#' @examples
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2, add_preds = TRUE)
#' post_pred_group(fit, grp = "Year", xlab = "Length (cm)")
#' @export
#'
post_pred_group <- function(mod, grp=NULL, xlab = 'Length (cm)'){

  pds <- mod$preds
  if (is.null(pds))
    stop("No predictive draws: refit with `fit_model(..., add_preds = TRUE)`")
  grps <- if (!is.null(grp)) syms(grp) else list()

  pps <- pds  %>%
    group_by(!!!grps, bin, .draw) %>%
    summarise(tot = sum(tot_by_bin), .pred = sum(.prediction), .groups = "drop_last") %>%
    group_by(!!!grps, .draw) %>%
    mutate( .pred = (.pred/sum(.pred))*sum(tot)) %>%
    group_by(!!!grps, bin, tot) %>%
    median_qi(.pred)

  p <- ggplot(pps) +
    geom_pointrange(aes(x=bin,y=.pred,ymin=.lower,ymax=.upper)) +
    geom_point(aes(x=bin,y=tot), col='skyblue') +
    xlab(xlab) +
    ylab('Prediction') +
    theme_compom() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))

  if(!is.null(grp)) p <- p+facet_wrap(vars(!!!grps),scales = 'free_y')

  return(p)
}

#' Plot bin-specific effects as multipliers
#'
#' Predicts the composition for each level of \code{grp}, holding all other
#' covariates at their mean (numeric) or mode (factor), and with the offset
#' fixed. Within each level and draw, predictions are converted to
#' proportions. Each proportion is then divided by the geometric mean of that
#' bin's proportions across levels of \code{grp} and draws. A multiplier above 1
#' indicates that the bin is relatively more frequent in that level than on
#' average.
#'
#' @param mod An object from fit_model
#' @param grp Character vector of the variable(s) to plot effects for (e.g.
#'   \code{"Year"}); must include \code{cvar} if one is used. Facets are by the
#'   first variable.
#' @param form Unused; retained for backwards compatibility.
#' @param grid Optional character vector of two variables in \code{grp} for a
#'   \code{facet_grid()} layout (rows, columns) instead of \code{facet_wrap()}.
#' @param cvar Name of the continuous covariate, if the model has one. It is
#'   held at its mean.
#' @param xlab x-axis label (the bin variable)
#' @return A ggplot of median multipliers and intervals per bin, faceted by
#'   the first variable in \code{grp} (or by \code{grid}).
#'
#' @importFrom tidybayes add_linpred_draws median_qi
#' @importFrom stats predict
#' @examples
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2)
#' Fx_plot(fit, grp = "Year")
#' @export
#'
Fx_plot <- function(mod, grp='Year', form='~(1|bin:Year)', grid=NULL, cvar=NULL,
                    xlab = 'Length bin (cm)'){

  if(!is.null(cvar)) grps <- syms(grp[(!grp %in% cvar)]) else grps <- syms(grp)

  #int <- as_draws_df(mod$mod, "b_Intercept") %>% dplyr::select(-.chain,-.iteration)
  helpers = any(grp %in% mod$data$bin)
  helper_vars = grp[grp %in% mod$data$bin]

  # Columns to exclude from factor conversion (coordinates for spatial models)
  exclude_cols <- paste0('^', grp, '$')

  if(any(grepl('TMB', attr(mod, "call")))) {
   #if(!is.null(coords)) exclude_col <- paste(exclude_col, paste0('^', coords, '$', collapse='|'), sep='|')
   exclude_cols <- paste(exclude_cols, paste0('^', mod$time, '$'), sep='|')
  }

  mean_or_mode <- function(x) {
    if(is.numeric(x)) mean(x) else factor(names(sort(-table(x)))[1])
  }

  #browser()
  preda_data <-
    mod$model$data %>%
    ungroup() %>%
    dplyr::select(-n, -tot_by_bin) %>%
    group_by(bin) %>%
    mutate(across(!matches(exclude_cols), ~mean_or_mode(.x))) %>%
    distinct() %>%
    mutate(n=100) %>% ungroup() %>%
    {if(!is.null(cvar)) mutate(., !!cvar:=mean(mod$data[[cvar]])) else .} %>%
    complete(nesting(!!!syms(grp),bin),fill=list(n=100))

  if(inherits(mod$model, 'brmsfit')){

    # brms re-evaluates factor(var) in the formula on newdata, which fails when a
    # factor held at its mode has a single level (e.g. Year when grp = "Season").
    # Pad with one row per fitted level of such factors; padded rows are
    # dropped after prediction (linear predictors are computed row-wise).
    fit_data <- mod$model$data
    collapsed <- names(preda_data)[vapply(names(preda_data), function(v)
      v %in% names(fit_data) && is.factor(fit_data[[v]]) &&
        !v %in% c(grp, "bin") &&
        dplyr::n_distinct(preda_data[[v]]) < nlevels(fit_data[[v]]),
      logical(1))]
    pad <- dplyr::bind_rows(lapply(collapsed, function(v) {
      lv <- setdiff(levels(fit_data[[v]]), as.character(preda_data[[v]]))
      p  <- preda_data[rep(1, length(lv)), , drop = FALSE]
      p[[v]] <- factor(lv, levels = levels(fit_data[[v]]))
      p
    }))
    preda_pad <- dplyr::bind_rows(
      dplyr::mutate(preda_data, .pad = FALSE),
      if (nrow(pad) > 0) dplyr::mutate(pad, .pad = TRUE)
    )

    preda <- preda_pad %>%
      add_linpred_draws(mod$model, allow_new_levels=T) %>%
      ungroup() %>%
      dplyr::filter(!.pad) %>%
      dplyr::select(-.pad)

  } else if(inherits(mod$model, 'sdmTMB')){

    # For spatial models, add coordinate columns using mean coordinates
    if(mod$spatial_varying && !is.null(mod$coords)){
      for(coord in mod$coords){
        if(!coord %in% names(preda_data) && coord %in% names(mod$data)){
          preda_data[[coord]] <- mean(mod$data[[coord]], na.rm=TRUE)
        }
      }
    }

    preda <- predict(mod$model,
                     newdata = preda_data,
                     offset = preda_data$n,
                     type = c("link"),
                     nsim = mod$nsim) %>%
      as.data.frame() %>%
      cbind(preda_data) %>%
      mutate(.row = 1:n()) %>%
      pivot_longer(cols=matches("V[0-9]+"), names_to = "iter", values_to = ".linpred") %>%
      mutate(.draw = as.numeric(gsub("V","",iter))) %>%
      ungroup()

  } else {stop("Model must be either brmsfit or sdmTMB")}

  # this is needed for when there are helper variables to run the splines over categories. We only want the combos of helpers
  # if(helpers){
  #
  #   preds <- c()
  #   for (v in helper_vars)  preds <- bind_rows(preds, preda[preda$bin==v, ] %>% filter(!!sym(v) == 1))
  #
  #   preds <- bind_rows(preds, preda[!preda$bin %in% helper_vars, ] %>% filter(across(helper_vars, ~. != 1)))
  #   preda <- preds
  # }
  #browser()
  preda <- preda %>%
    group_by(!!!syms(grps),.draw) %>%
    mutate(lp = exp(.linpred)/sum(exp(.linpred))) %>%
    group_by(bin) %>%
    mutate(ml=gmean(lp)) %>%
    group_by(bin,!!!syms(grps)) %>%
    median_qi(pred = lp/ml)


  p <- ggplot(preda)  + geom_hline(yintercept = 1,linetype=2,alpha=0.5) +
    geom_pointrange(aes(x=bin,y=pred,ymin=.lower,ymax=.upper)) +
    xlab(xlab) +
    ylab('Multiplier') +
    theme_compom() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))
  #browser()
  if(is.null(grid))  p <- p+facet_wrap(facets = vars(!!grps[[1]]),scales = 'free_y')
  if(!is.null(grid)) p <- p+facet_grid(rows = vars(!!as.symbol(grid[1])),cols = vars(!!as.symbol(grid[2])),scales = 'free_y')

  return(p)
}

#' Scale predicted compositions to catch or to an index population
#'
#' Predicts numbers-at-bin from a fitted model and weights the predictions.
#' The result is a set of draws of numbers-at-bin per stratum (or grid cell),
#' which can be summed over any subset of \code{grps} and normalised to obtain
#' compositions (e.g. with \code{scaled_comp_plot()}).
#'
#' \code{mode} selects between two uses.
#' \describe{
#'   \item{\code{"catch"} (default)}{The composition of the catch. Each stratum
#'     in \code{scale_df} is predicted at its own level of every model term,
#'     including the nuisance random effects (\code{re_form}) and continuous
#'     covariates, and the prediction is weighted by the catch of the stratum
#'     (\code{predvar}). Every fixed-effect and nuisance variable must
#'     therefore be in \code{grps}.}
#'   \item{\code{"index"}}{The composition of the population that a CPUE index
#'     describes, predicted as a CPUE analysis would predict the index.
#'     Catchability terms (\code{catchability}) are set to a reference. All
#'     other terms are predicted at the values in \code{grid}. For non-spatial
#'     models, \code{grid} has one row per index time step (e.g. Year). For
#'     spatial models, \code{grid} has one row per cell and time step, with the
#'     coordinates, the time column and the density covariates. Cells are
#'     weighted by \code{weight} (equal weights if \code{NULL}), which can come
#'     from a separate CPUE or abundance model.}
#' }
#'
#' In index mode, every fixed-effect variable and continuous covariate of the
#' model is either listed in \code{catchability} or is a column of \code{grid}.
#' The function stops otherwise, so that the role of each term is chosen by the
#' user. Random effects (\code{re_form}) are catchability terms unless they are
#' listed in \code{density}, in which case they are evaluated at the values in
#' \code{grid}. With the TMB backend, either all or none of the random effects
#' can be kept.
#'
#' Reference values for catchability terms are the most frequent level (factors)
#' or the mean (continuous covariates) in the fitted data, unless given in
#' \code{reference}. Random effects are set to the population level.
#'
#' In both modes:
#' \itemize{
#'   \item Continuous covariates (\code{cvars}) are evaluated at their mean
#'     within each stratum or cell of the data frame supplied.
#'   \item Spatial and spatio-temporal fields (spatial sdmTMB fits) are
#'     evaluated at the coordinates in the data frame supplied.
#' }
#'
#' Nuisance levels that were not in the fitted data (e.g. unsampled vessels)
#' are handled differently by the backends. brms uses, in each draw, the
#' posterior draw of a randomly chosen fitted level
#' (\code{sample_new_levels = "uncertainty"}). sdmTMB predicts such rows at
#' the population level of all nuisance random effects, and a warning is
#' issued.
#'
#' @param scale_df data frame used for catch weighting, with one or more rows
#'   per stratum, containing \code{grps}, \code{predvar}, any \code{cvars}, and
#'   (spatial fits) the coordinate and time columns.
#' @param predvar Column in \code{scale_df} used for scaling (e.g. catch)
#' @param grps Grouping variables defining the strata (catch mode) or the
#'   output time steps and other prediction variables (index mode). Must
#'   include the \code{ffx_form} variables (except catchability terms) and, for
#'   catch weighting, the \code{re_form} variables.
#' @param fit A \code{ComPoM} object from \code{fit_model()}
#' @param iters Number of posterior draws (brms only; default uses all stored
#'   draws). For sdmTMB fits the number of draws is set by \code{fit_model()}.
#' @param pgrid Not used; reserved for future use.
#' @param remove_nuisance Deprecated. \code{FALSE} is \code{mode = "catch"}.
#'   \code{TRUE} is \code{mode = "index"} with the random effects as
#'   catchability terms, \code{scale_df} as \code{grid} and \code{predvar} as
#'   \code{weight}.
#' @param mode \code{"catch"} or \code{"index"}; see Details.
#' @param catchability Character vector of model terms that are catchability
#'   covariates, set to a reference (index mode).
#' @param density Optional character vector of terms declared as density terms
#'   (index mode). They must be columns of \code{grid}. Random effects listed
#'   here are kept.
#' @param reference Named list of reference values for catchability terms, e.g.
#'   \code{list(Gear = "trawl", sst = 15)}.
#' @param grid Prediction grid for \code{mode = "index"}; see Details. Defaults
#'   to \code{scale_df}.
#' @param weight Column of \code{grid} with the weight of each row (index
#'   mode). \code{NULL} gives equal weights.
#' @return A data frame of draws (\code{.draw}) with predicted numbers-at-bin
#'   (\code{tot_by_bin}) per stratum or cell and bin. The weight per row (sum
#'   of the weight column) is in \code{n}.
#'
#' @importFrom tidybayes add_predicted_draws
#' @importFrom stats predict
#' @seealso \code{scaled_comp_plot()}, \code{scaled_ridge_plot()}
#' @examples
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2)
#'
#' # catch weighting: composition of the total catch of the fleet
#' sc <- scale_comps(ex_catch, predvar = "catch", fit = fit,
#'                   grps = c("Year", "Season", "Vessel"))
#' head(sc)
#' scaled_comp_plot(sc, grps = "Year")
#'
#' # index population: Season and Vessel are catchability terms; one cell per year
#' sc_pop <- scale_comps(fit = fit, mode = "index", grps = "Year",
#'                       grid = data.frame(Year = factor(2015:2020)))
#' scaled_comp_plot(sc_pop, grps = "Year")
#' @export
#'
scale_comps <- function(scale_df,
                        predvar    = 'catch',
                        fit        = NULL,
                        grps       = NULL,
                        iters      = NULL,
                        pgrid      = NULL,
                        remove_nuisance = FALSE,
                        mode       = c("catch", "index"),
                        catchability = NULL,
                        density    = NULL,
                        reference  = NULL,
                        grid       = NULL,
                        weight     = NULL) {

  if (is.null(fit)) stop("'fit' must be a ComPoM object from fit_model()")
  mode_given <- !missing(mode)
  mode <- match.arg(mode)

  # Deprecated interface: remove_nuisance = TRUE
  if (isTRUE(remove_nuisance)) {
    if (mode_given && mode == "catch")
      stop("remove_nuisance = TRUE conflicts with mode = 'catch'")
    rlang::warn("`remove_nuisance` is deprecated; use `mode = \"index\"`.",
                .frequency = "once", .frequency_id = "ComPoM_remove_nuisance")
    mode   <- "index"
    weight <- if (is.null(weight)) predvar else weight
  }

  if (mode == "index") {
    if (is.null(grid)) {
      if (missing(scale_df)) stop("`grid` must be supplied for mode = 'index'")
      grid <- scale_df
    }
    grid <- as.data.frame(grid)
    if (!is.null(weight) && !weight %in% names(grid))
      stop("`weight` column '", weight, "' is not in `grid`")
    grid$.w <- if (is.null(weight)) 1 else grid[[weight]]
    scale_df <- grid
    predvar  <- ".w"
  }

  bin_levels    <- unique(fit$data$bin)
  bin_is_factor <- is.factor(fit$data$bin)
  cvars         <- fit$cvars  # continuous spline covariates (may be NULL)

  # ── helper: build newdata with one row per group x bin ──────────────────────
  # If cvars are present they must be columns in scale_df; mean per group is
  # used so predictions represent average oceanographic conditions per stratum.
  build_newdata <- function(df, group_vars) {
    nd <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        n = sum(!!rlang::sym(predvar), na.rm = TRUE),
        dplyr::across(dplyr::any_of(cvars), mean),
        .groups = "drop"
      ) |>
      dplyr::mutate(bin = if (bin_is_factor) factor(bin_levels[1], levels = bin_levels)
                         else bin_levels[1]) |>
      tidyr::complete(tidyr::nesting(!!!rlang::syms(group_vars)), bin,
                      fill = list(n = 0)) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::mutate(
        n   = sum(n),
        bin = if (bin_is_factor) bin else as.numeric(as.character(bin))
      ) |>
      dplyr::filter(n > 0) |>
      dplyr::ungroup()

    # Re-join cvars after complete() (complete() fills them with NA)
    if (!is.null(cvars)) {
      cvar_means <- df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
        dplyr::summarise(dplyr::across(dplyr::any_of(cvars), mean), .groups = "drop")
      nd <- nd |>
        dplyr::select(-dplyr::any_of(cvars)) |>
        dplyr::left_join(cvar_means, by = group_vars)
      # cbin = numeric version of bin, required by t2(sst, cbin) smoother
      nd$cbin <- as.numeric(as.character(nd$bin))
    }
    nd
  }

  # ── helper: factor-coerce newdata for TMB (preserve coords, time, bin, cvars) ──
  coerce_tmb <- function(df) {
    excl <- ".*_by_.*|^n$|^cbin$"
    if (!bin_is_factor) excl <- paste(excl, "^bin$", sep = "|")
    if (!is.null(fit$coords))
      excl <- paste(excl, paste0("^", fit$coords, "$", collapse = "|"), sep = "|")
    if (!is.null(fit$time))
      excl <- paste(excl, paste0("^", fit$time, "$"), sep = "|")
    if (!is.null(cvars))
      excl <- paste(excl, paste0("^", cvars, "$", collapse = "|"), sep = "|")
    df |> dplyr::mutate(dplyr::across(-dplyr::matches(excl), as.factor))
  }

  # ── helper: pivot TMB simulation matrix to long draws ───────────────────────
  # sim_mat contains predicted counts at offset=0 (i.e. proportions * 1).
  # Multiply by catch (nd$n) per row to get catch-scaled compositions.
  tmb_to_long <- function(sim_mat, newdata) {
    # predicted values are proportional (offset=log(1)=0); scale by catch
    sim_mat <- sim_mat * newdata$n
    sim_mat |>
      as.data.frame() |>
      cbind(newdata) |>
      dplyr::mutate(.row = dplyr::row_number()) |>
      tidyr::pivot_longer(
        cols      = dplyr::matches("V[0-9]+"),
        names_to  = "iter",
        values_to = "tot_by_bin"
      ) |>
      dplyr::mutate(.draw = as.numeric(gsub("V", "", iter))) |>
      dplyr::ungroup()
  }

  group_vars <- if (fit$spatial_varying) unique(c(grps, fit$coords)) else grps
  nd <- build_newdata(scale_df, group_vars)

  # ── term roles ──────────────────────────────────────────────────────────────
  # Catch mode: every term is evaluated at the stratum's own level, so every
  # nuisance variable must be in grps (and scale_df).
  # Index mode: catchability terms are set to a reference; all other terms are
  # predicted at the values in the grid. Random effects are catchability terms
  # unless listed in `density`.
  re_vars <- fit$re_vars
  keep_re_vars <- re_vars
  if (mode == "catch") {
    missing_re <- setdiff(re_vars, names(nd))
    if (length(missing_re) > 0)
      stop("Catch weighting (mode = 'catch') predicts each stratum at its own ",
           "level of every nuisance variable; add ",
           paste0("'", missing_re, "'", collapse = ", "),
           " to `grps` and `scale_df`, or use mode = 'index'.")
  } else {
    model_terms <- unique(c(fit$ffx_vars, cvars))
    bad_cat <- setdiff(catchability, c(model_terms, re_vars))
    if (length(bad_cat) > 0)
      stop("`catchability` terms not in the model: ", paste(bad_cat, collapse = ", "))
    in_grps <- intersect(catchability, group_vars)
    if (length(in_grps) > 0)
      stop("Catchability terms are set to a reference and must not be in `grps`: ",
           paste(in_grps, collapse = ", "))
    if (length(intersect(catchability, density)) > 0)
      stop("Terms cannot be both catchability and density: ",
           paste(intersect(catchability, density), collapse = ", "))
    bad_den <- setdiff(density, names(nd))
    if (length(bad_den) > 0)
      stop("`density` terms must be columns of `grid`: ", paste(bad_den, collapse = ", "))
    unclassified <- setdiff(model_terms, c(catchability, names(nd)))
    if (length(unclassified) > 0)
      stop("Model terms with no role in index mode: ", paste(unclassified, collapse = ", "),
           ". List them in `catchability`, or add them to `grid`.")
    keep_re_vars <- intersect(re_vars, density)

    # reference values for catchability fixed effects and covariates
    for (v in intersect(catchability, model_terms)) {
      x   <- fit$data[[v]]
      val <- if (!is.null(reference[[v]])) reference[[v]] else
        if (is.numeric(x)) mean(x, na.rm = TRUE) else names(sort(-table(x)))[1]
      nd[[v]] <- if (is.numeric(x)) as.numeric(val) else
        factor(as.character(val), levels = levels(x))
    }
  }
  if (length(keep_re_vars) > 0 && length(keep_re_vars) < length(re_vars) &&
      inherits(fit$model, "sdmTMB"))
    stop("With the TMB backend, either all or none of the random effects can be kept: ",
         "list all of ", paste(re_vars, collapse = ", "), " in `density`, or none.")
  keep_re <- length(keep_re_vars) > 0
  re_formula <- if (!keep_re) NA else if (length(keep_re_vars) == length(re_vars)) NULL else
    stats::as.formula(paste("~", paste0("(1|bin:", keep_re_vars, ")", collapse = " + ")))

  if (inherits(fit$model, "brmsfit")) {

    # Levels not seen in the fit (e.g. unsampled vessels) take, in each draw,
    # the posterior draw of a randomly chosen fitted level
    # (brms sample_new_levels = "uncertainty").
    return(tidybayes::add_predicted_draws(
      nd, fit$model,
      allow_new_levels = TRUE,
      ndraws           = iters,
      value            = "tot_by_bin",
      re_formula       = re_formula
    ))
  }

  if (!inherits(fit$model, "sdmTMB"))
    stop("Model must be either brmsfit or sdmTMB")

  # predict_data carries the model columns; nd keeps the stratum values for output
  nd <- coerce_tmb(nd)
  predict_data <- nd
  for (v in re_vars) {
    lev <- if (keep_re) as.character(nd[[v]]) else levels(fit$data[[v]])[1]
    if (keep_re) {
      new_lev <- setdiff(unique(lev), levels(fit$data[[v]]))
      if (length(new_lev) > 0)
        warning("Levels of '", v, "' not in the fitted data: ",
                paste(utils::head(new_lev, 5), collapse = ", "),
                if (length(new_lev) > 5) ", ..." else "",
                ". sdmTMB predicts these rows at the population level of ALL nuisance ",
                "random effects (re_form_iid = NA for those rows).", call. = FALSE)
    }
    predict_data[[v]]                 <- factor(lev)
    predict_data[[paste0("bin:", v)]] <- factor(paste(nd$bin, lev, sep = "_"))
  }

  if (fit$spatial_varying) {
    predict_data <- sdmTMB::make_category_svc(
      data                    = predict_data,
      category_column         = "bin",
      time_column             = fit$time,
      share_spatial_sd        = TRUE,
      share_spatiotemporal_sd = TRUE
    )$data_expanded
  }

  sim_mat <- predict(
    fit$model,
    type             = "response",
    newdata          = predict_data,
    offset           = rep(0, nrow(predict_data)),
    nsim             = fit$nsim,
    re_form_iid      = if (keep_re) NULL else NA,
    allow_new_levels = TRUE
  )

  tmb_to_long(sim_mat, nd)
}


#' Plot scaled compositions
#'
#' Sums the draws from \code{scale_comps()} over all strata within levels of
#' \code{grps} and plots the resulting compositions (as proportions or as
#' catch-at-bin) with uncertainty intervals. Optional reference compositions
#' (e.g. raw samples or a known truth) can be overlaid.
#'
#' @param scaled_comp output of \code{scale_comps()}
#' @param grps grouping variables to aggregate to (e.g. \code{"Year"})
#' @param scaled If \code{TRUE}, plot proportions; if \code{FALSE}, plot
#'   catch-at-bin.
#' @param scales Passed to \code{facet_wrap()}/\code{facet_grid()}. Panels are
#'   a \code{facet_wrap()} by \code{grps} when it has one element, and a
#'   \code{facet_grid()} with the first element in rows and the remaining
#'   elements in columns otherwise.
#' @param comp_are,comp_are2 Optional data frames with \code{grps}, \code{bin}
#'   and \code{cvar}, overlaid as blue and orange points respectively. Counts
#'   are converted to proportions within \code{grps} when \code{scaled = TRUE}.
#' @param cvar Column in \code{comp_are}/\code{comp_are2} holding the counts
#'   or catch per bin.
#' @param xlab x-axis label (the bin variable)
#' @param comp_label,comp2_label Legend labels for the points from
#'   \code{comp_are} (skyblue) and \code{comp_are2} (orange).
#' @return A ggplot with the mean and 95\% quantile interval of the
#'   composition (or catch-at-bin) per bin.
#'
#' @examples
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2)
#' sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
#'
#' scaled_comp_plot(sc, grps = "Year")
#'
#' # overlay the raw sample composition
#' raw <- dplyr::rename(ex_comp, bin = length)
#' scaled_comp_plot(sc, grps = "Year", comp_are = raw, cvar = "count")
#' @export
#'
scaled_comp_plot <- function(scaled_comp=NULL,
                             grps = NULL,
                             scaled=T,
                             scales = "free_y",
                             comp_are = NULL,
                             comp_are2 = NULL,
                             cvar='tcatch',
                             xlab = 'Length bin (cm)',
                             comp_label = 'Reference',
                             comp2_label = 'Reference 2'){

  if(scaled == T) {

    preda <- scaled_comp %>%
      group_by(across(all_of(grps) ), bin, .draw) %>%
      summarise(tot = sum(tot_by_bin), .groups = "drop_last") %>%
      group_by(across(all_of(grps) ), .draw) %>%
      mutate(prop = tot/sum(tot)) %>%
      group_by(across(all_of(grps) ), bin) %>%
      mean_qi(prop) %>%
      filter(!is.na(prop)) %>%
      mutate(bin = as.numeric(as.character(bin)))

    if(!is.null(comp_are)){
      comp <- comp_are %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(tot = sum(!!sym(cvar)), .groups = "drop") %>%
        group_by(across(all_of(grps) )) %>%
        mutate(prop = tot/sum(tot)) %>%
        filter(!is.na(prop), prop > 0) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

    if(!is.null(comp_are2)){
      comp2 <- comp_are2 %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(tot = sum(!!sym(cvar)), .groups = "drop") %>%
        group_by(across(all_of(grps) )) %>%
        mutate(prop = tot/sum(tot)) %>%
        filter(!is.na(prop), prop > 0) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

  } else {

    preda <- scaled_comp %>%
      group_by(across(all_of(grps) ), bin, .draw) %>%
      summarise(tot = sum(tot_by_bin)) %>%
      group_by(across(all_of(grps) ), bin) %>%
      mean_qi(prop = tot)  %>%
      filter(!is.na(prop)) %>%
      mutate(bin = as.numeric(as.character(bin)))

    if(!is.null(comp_are)){

      comp <- comp_are %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(prop = sum(!!sym(cvar))) %>%
        filter(!is.na(prop),sum(prop)>100) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }

    if(!is.null(comp_are2)){

      comp2 <- comp_are2 %>%
        group_by(across(all_of(grps) ), bin) %>%
        summarise(prop = sum(!!sym(cvar))) %>%
        filter(!is.na(prop),sum(prop)>100) %>%
        mutate(bin = as.numeric(as.character(bin)))
    }
  }


  p <- ggplot( preda ) +
    geom_pointrange(aes(x=bin,y=prop,ymin=.lower,ymax=.upper)) +
    {if(!is.null(comp_are)) geom_point(aes(x=bin,y=prop,colour=comp_label), data=comp)} +
    {if(!is.null(comp_are2)) geom_point(aes(x=bin,y=prop,colour=comp2_label), data=comp2)} +
    scale_colour_manual(
      values = c(setNames('skyblue', comp_label)[!is.null(comp_are)],
                 setNames('orange', comp2_label)[!is.null(comp_are2)]),
      breaks = c(comp_label, comp2_label)[c(!is.null(comp_are), !is.null(comp_are2))],
      name = NULL) +
    scale_x_continuous(breaks = int_breaks) +
    xlab(xlab) +
    ylab('Prediction') +
    theme_compom() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8),
          legend.position = 'bottom')

  if(length(grps)==1) p <- p+facet_wrap(facets = vars(!!as.symbol(grps)),scales = scales,drop = T,ncol = 4)
  if(length(grps)>=2) p <- p+facet_grid(rows = vars(!!as.symbol(grps[1])),cols = vars(!!!syms(grps[2:length(grps)])),scales = scales)

  p
}


#' Plot scaled compositions as ridgeline plot
#'
#' Sums the draws from \code{scale_comps()} over strata within levels of
#' \code{grps}, converts them to proportions and plots the median composition
#' of each level as a ridgeline, stacked along \code{plot_y}.
#'
#' @param scaled_comp Data frame from \code{scale_comps()}.
#' @param grps Character vector of grouping variables to aggregate to; must
#'   contain \code{plot_y} (and \code{plot_facet}, if used).
#' @param plot_y Variable name (string) to stack on the y-axis (e.g.
#'   \code{"Year"}).
#' @param plot_facet Optional variable name (string) to facet by. \code{NULL} =
#'   no facet.
#' @param scales Passed to \code{facet_wrap()}, default \code{"free_y"}.
#' @param bin_range Optional numeric vector \code{c(min, max)} to restrict
#'   plotted bins.
#' @param xlab x-axis label (the bin variable)
#' @return A ggplot.
#'
#' @import cowplot
#' @importFrom ggridges geom_ridgeline
#' @examples
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2)
#' sc <- scale_comps(ex_catch, "catch", fit, grps = c("Year", "Season", "Vessel"))
#' scaled_ridge_plot(sc, grps = "Year", plot_y = "Year")
#' @export
#'
scaled_ridge_plot <- function(scaled_comp = NULL,
                              grps        = NULL,
                              plot_y      = NULL,
                              plot_facet  = NULL,
                              scales      = "free_y",
                              bin_range   = NULL,
                              xlab        = "Length bin (cm)") {

  if (is.null(plot_y)) stop("`plot_y` must be supplied as a column name string")

  sdat <- scaled_comp |>
    dplyr::ungroup() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grps)), .draw, bin) |>
    dplyr::summarise(n = sum(tot_by_bin), .groups = "drop") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grps)), .draw) |>
    dplyr::mutate(p = n / sum(n)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grps)), bin) |>
    tidybayes::median_qi(p) |>
    dplyr::mutate(bin = as.numeric(as.character(bin)))

  if (!is.null(bin_range)) {
    sdat <- dplyr::filter(sdat, bin >= bin_range[1], bin <= bin_range[2])
  }

  p <- ggplot2::ggplot(
    sdat,
    ggplot2::aes(
      x      = bin,
      y      = !!rlang::sym(plot_y),
      height = p
    )
  ) +
    ggridges::geom_ridgeline(
      scale     = 8,
      fill      = "skyblue",
      alpha     = 0.6,
      color     = "grey30",
      linewidth = 0.3
    ) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(plot_y) +
    theme_compom()

  if (!is.null(plot_facet)) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym(plot_facet)),
      scales = scales
    )
  }

  return(p)
}




#' Plot marginal smooth effects of a continuous covariate by bin
#'
#' Plots the centred smooth effect of a continuous covariate with 95\%
#' confidence intervals for each bin level. Works for both sdmTMB and brms
#' backends, and for both \code{spl = "by_factor"} (\code{s(cvar, by=bin)})
#' and \code{spl = "2D"} (\code{t2(cvar, cbin)}) spline types.
#' Fixed effects and IID random effects are marginalised out so only the
#' covariate smooth contribution is shown.
#'
#' @param fit A \code{ComPoM} object from \code{fit_model()} with
#'   \code{cvars} set.
#' @param cvar Character. Name of the continuous covariate to plot. Defaults
#'   to the first element of \code{fit$cvars}.
#' @param ref_ffx Named list of reference values for the fixed-effect factor
#'   variables (e.g. \code{list(yy = 2015)}). Defaults to the first observed
#'   value of each fixed-effect variable detected from the model formula.
#' @param ref_re Character. Reference level of the nuisance RE variable.
#'   Defaults to the first level of the first RE variable.
#' @param n_grid Integer. Number of covariate values in the prediction grid.
#'   Default 50.
#' @return A ggplot of the centred smooth effect (log scale) per bin.
#' @examples
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip", "sst"), "length", "count")
#' dat$sst <- as.numeric(as.character(dat$sst))
#' fit <- fit_model("Year", "Season", cvars = "sst", spl = "by_factor", knots = "4",
#'                  data = dat, backend = "TMB", iter = 300, warmup = 100, thin = 2)
#' cvar_smooth_plot(fit, "sst")
#' @export
cvar_smooth_plot <- function(fit,
                             cvar    = NULL,
                             ref_ffx = NULL,
                             ref_re  = NULL,
                             n_grid  = 50L) {

  if (is.null(cvar)) cvar <- fit$cvars[1]
  if (is.null(cvar)) stop("No continuous covariate found in fit$cvars")

  is_brms  <- inherits(fit$model, "brmsfit")
  mod_data <- if (is_brms) fit$model$data else fit$model$data
  bins     <- levels(fit$data$bin)

  # ── Detect fixed-effect factor variable names from factor(var):bin columns ──
  ffy_cols <- grep("^factor\\(", names(mod_data), value = TRUE)
  ffx_vars <- unique(gsub("factor\\(([^)]+)\\).*", "\\1", ffy_cols))

  # ── Default reference: first observed value of each ffx variable ─────────
  if (is.null(ref_ffx)) {
    ref_ffx <- lapply(ffx_vars, function(v) {
      if (v %in% names(mod_data)) sort(unique(mod_data[[v]]))[1] else NULL
    })
    names(ref_ffx) <- ffx_vars
  }

  # ── Default RE reference level ────────────────────────────────────────────
  if (is.null(ref_re) && !is.null(fit$re_vars))
    ref_re <- levels(fit$data[[fit$re_vars[1]]])[1]

  cvar_seq <- seq(
    min(mod_data[[cvar]], na.rm = TRUE),
    max(mod_data[[cvar]], na.rm = TRUE),
    length.out = n_grid
  )

  # ── Build prediction grid ─────────────────────────────────────────────────
  # Template row per bin: find a row matching the reference ffx values
  grid <- purrr::map_dfr(bins, function(b) {
    row <- mod_data
    for (v in names(ref_ffx))
      if (!is.null(ref_ffx[[v]]) && v %in% names(row))
        row <- row[row[[v]] == ref_ffx[[v]], , drop = FALSE]
    row <- row[row$bin == b, , drop = FALSE]
    if (nrow(row) == 0) row <- mod_data[mod_data$bin == b, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    row <- row[1, , drop = FALSE]

    purrr::map_dfr(cvar_seq, function(s) {
      r         <- row
      r[[cvar]] <- s
      # cbin is the numeric version of bin used by t2(cvar, cbin)
      if ("cbin" %in% names(r))
        r$cbin <- as.numeric(as.character(r$bin[1]))
      r
    })
  })

  # Zero all factor(var):bin dummy columns — marginalise fixed effects
  present_ffy <- ffy_cols[ffy_cols %in% names(grid)]
  grid[present_ffy] <- 0L

  # ── Set RE columns to reference level ─────────────────────────────────────
  if (!is.null(fit$re_vars) && !is.null(ref_re)) {
    for (v in fit$re_vars) {
      grid[[v]] <- factor(ref_re, levels = levels(fit$data[[v]]))
      col       <- paste0("bin:", v)
      if (col %in% names(mod_data))
        grid[[col]] <- factor(paste(grid$bin, ref_re, sep = "_"),
                              levels = levels(mod_data[[col]]))
    }
  }

  # ── Predict — backend-specific ────────────────────────────────────────────
  if (is_brms) {

    # posterior_linpred returns a draws x rows matrix (linear scale, offset=0)
    grid$n <- 1L
    draws <- brms::posterior_linpred(
      fit$model,
      newdata          = grid,
      re_formula       = NA,
      allow_new_levels = TRUE
    )  # draws x (n_bins * n_grid) matrix

    # Centre each draw within bin, then exponentiate → multiplier per draw
    n_bins  <- length(bins)
    n_grid_pts <- nrow(grid) / n_bins   # points per bin

    plot_df <- purrr::map_dfr(seq_along(bins), function(bi) {
      # rows are bin-major: all n_grid_pts rows for bin bi are consecutive
      idx    <- ((bi - 1) * n_grid_pts + 1):(bi * n_grid_pts)
      d_bin  <- draws[, idx, drop = FALSE]         # draws x n_grid_pts
      # Centre per draw: subtract per-draw mean so multiplier = 1 at mean SST
      d_c    <- exp(d_bin - rowMeans(d_bin))
      tibble::tibble(
        bin    = bins[bi],
        !!cvar := cvar_seq,
        est_c  = colMeans(d_c),
        lo     = apply(d_c, 2, quantile, 0.025),
        hi     = apply(d_c, 2, quantile, 0.975)
      )
    }) |>
      dplyr::mutate(bin = factor(bin, levels = bins))

  } else {

    # sdmTMB: draw nsim samples from the joint precision, offset = 0
    sim_mat <- predict(
      fit$model,
      newdata     = grid,
      offset      = rep(0, nrow(grid)),
      re_form     = NA,
      re_form_iid = NA,
      nsim        = fit$nsim
    )  # rows x nsim matrix (linear scale)

    n_bins     <- length(bins)
    n_grid_pts <- nrow(grid) / n_bins

    plot_df <- purrr::map_dfr(seq_along(bins), function(bi) {
      # rows are bin-major: consecutive block for bin bi
      idx   <- ((bi - 1) * n_grid_pts + 1):(bi * n_grid_pts)
      d_bin <- t(sim_mat[idx, , drop = FALSE])   # nsim x n_grid_pts
      # Centre per draw, then exponentiate → multiplier on response scale
      d_c   <- exp(d_bin - rowMeans(d_bin))
      tibble::tibble(
        bin    = bins[bi],
        !!cvar := cvar_seq,
        est_c  = colMeans(d_c),
        lo     = apply(d_c, 2, quantile, 0.025),
        hi     = apply(d_c, 2, quantile, 0.975)
      )
    }) |>
      dplyr::mutate(bin = factor(bin, levels = bins))
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[cvar]], y = est_c)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                         fill = "skyblue", alpha = 0.25) +
    ggplot2::geom_line(colour = "skyblue", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::facet_wrap(~bin, labeller = ggplot2::label_both, scales = "free_y") +
    ggplot2::labs(
      x     = cvar,
      y     = "Multiplier (response scale, centred)"
    ) +
    theme_compom()
}


#' Plot predicted spatial field by bin
#'
#' Produces a faceted tile map showing the spatial random field contribution
#' per bin, centred by subtracting each bin mean so the map shows spatial
#' deviation from the average availability of that bin size class.
#'
#' @param fit A ComPoM object from fit_model() with a spatial mesh.
#' @param cvar Character. Name of the continuous covariate. If supplied the
#'   spline is evaluated at the mean observed value. Default NULL.
#' @param ref_ffx Named list of reference values for the fixed-effect factor
#'   variables (e.g. list(yy = 2015)). Defaults to the first observed value of
#'   each fixed-effect variable detected from the model formula.
#' @param ref_re Character. Reference level of the nuisance RE. Default is the first level.
#' @param lon_step Numeric. Longitude grid resolution. Default 5.
#' @param lat_step Numeric. Latitude grid resolution. Default 5.
#' @param ncol Integer. Number of facet columns. Default 4.
#' @return A ggplot with one tile map per bin.
#' @examples
#' \donttest{
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip", "x", "y"),
#'                  "length", "count")
#' dat$x <- as.numeric(as.character(dat$x))
#' dat$y <- as.numeric(as.character(dat$y))
#' mesh <- sdmTMB::make_mesh(dat, c("x", "y"), cutoff = 25)
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  mesh = mesh, time = "Year", coords = c("x", "y"),
#'                  iter = 300, warmup = 100, thin = 2)
#' spatial_field_plot(fit, lon_step = 10, lat_step = 10)
#' }
#' @export
spatial_field_plot <- function(fit,
                               cvar     = NULL,
                               ref_ffx  = NULL,
                               ref_re   = NULL,
                               lon_step = 5,
                               lat_step = 5,
                               ncol     = 4L) {

  if (!fit$spatial_varying)
    stop("`spatial_field_plot()` requires a spatial sdmTMB model")

  mod_data  <- fit$model$data
  bins      <- levels(fit$data$bin)
  lon_col   <- fit$coords[1]
  lat_col   <- fit$coords[2]

  # Detect fixed-effect factor variable names from factor(var):bin columns
  ffy_cols <- grep("^factor\\(", names(mod_data), value = TRUE)
  ffx_vars <- unique(gsub("factor\\(([^)]+)\\).*", "\\1", ffy_cols))

  # Default reference: first observed value of each ffx variable
  if (is.null(ref_ffx)) {
    ref_ffx <- lapply(ffx_vars, function(v) {
      if (v %in% names(mod_data)) sort(unique(mod_data[[v]]))[1] else NULL
    })
    names(ref_ffx) <- ffx_vars
  }

  # Default RE reference level
  if (is.null(ref_re) && !is.null(fit$re_vars))
    ref_re <- levels(fit$data[[fit$re_vars[1]]])[1]

  lon_seq   <- seq(min(mod_data[[lon_col]]), max(mod_data[[lon_col]]), by = lon_step)
  lat_seq   <- seq(min(mod_data[[lat_col]]), max(mod_data[[lat_col]]), by = lat_step)
  mean_cvar <- if (!is.null(cvar)) mean(mod_data[[cvar]], na.rm = TRUE) else NULL

  # Build spatial prediction grid
  grid <- tidyr::expand_grid(
    !!lon_col := lon_seq,
    !!lat_col := lat_seq,
    bin = factor(bins, levels = bins)
  ) |>
    dplyr::mutate(
      n          = 1L,
      tot_by_bin = 0L,
      cbin       = as.numeric(as.character(bin))
    )

  # Set fixed-effect factor variables to reference level
  for (v in names(ref_ffx)) {
    if (!is.null(ref_ffx[[v]])) {
      grid[[v]] <- if (v %in% names(fit$data) && is.factor(fit$data[[v]]))
        factor(ref_ffx[[v]], levels = levels(fit$data[[v]])) else ref_ffx[[v]]
    }
  }

  if (!is.null(cvar)) grid[[cvar]] <- mean_cvar

  # Add bin dummy columns
  for (b in bins)
    grid[[paste0("bin", b)]] <- as.integer(grid$bin == b)

  # Zero all factor(var):bin dummy columns to marginalise fixed effects
  for (col in ffy_cols) grid[[col]] <- 0L

  # Set RE columns to reference level
  if (!is.null(fit$re_vars) && !is.null(ref_re)) {
    for (v in fit$re_vars) {
      grid[[v]] <- factor(ref_re, levels = levels(fit$data[[v]]))
      col       <- paste0("bin:", v)
      if (col %in% names(mod_data))
        grid[[col]] <- factor(paste(grid$bin, ref_re, sep = "_"),
                              levels = levels(mod_data[[col]]))
    }
  }

  preds <- predict(
    fit$model,
    newdata     = grid,
    offset      = rep(0, nrow(grid)),
    re_form_iid = NA,
    se_fit      = FALSE
  )

  lim <- preds |>
    dplyr::group_by(bin) |>
    dplyr::mutate(est_c = est - mean(est, na.rm = TRUE)) |>
    dplyr::pull(est_c) |> abs() |> max(na.rm = TRUE)

  plot_df <- preds |>
    dplyr::group_by(bin) |>
    dplyr::mutate(est_c = est - mean(est, na.rm = TRUE)) |>
    dplyr::ungroup()

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(x = .data[[lon_col]],
                               y = .data[[lat_col]],
                               fill = est_c)) +
    ggplot2::geom_tile(width = lon_step, height = lat_step) +
    ggplot2::facet_wrap(~bin, labeller = ggplot2::label_both, ncol = ncol) +
    ggplot2::scale_fill_viridis_c(
      option = "cividis", limits = c(-lim, lim),
      name    = "Spatial\ndeviation\n(log scale)"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x     = lon_col,
      y     = lat_col
    ) +
    theme_compom() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

# Axis breaks at whole numbers (bins are lengths or ages, not fractions)
int_breaks <- function(lims) unique(round(pretty(lims)))

