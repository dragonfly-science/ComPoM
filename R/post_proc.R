#' Posterior predictive check by group
#'
#' Compares observed counts per bin with posterior predictive (brms) or
#' simulated (sdmTMB) counts, summed within levels of \code{grp}. Predicted
#' totals are rescaled to the observed total within each group and draw, so
#' the check is on the composition rather than on the total.
#' Requires a fit from \code{fit_model(..., add_preds = TRUE)}.
#'
#' @param mod An object from fit_model
#' @param grp group to do PPC over
#' @param xlab label for composition plot
#' @return A ggplot with predicted medians and intervals (black) and observed
#'   counts (blue), faceted by \code{grp}.
#'
#' @import cowplot
#' @import ggplot2
#' @export
#'
post_pred_group <- function(mod, grp=NULL, xlab = 'Length (cm)'){

  pds <- mod$preds
  if(!is.null(grp)) grps <- sym(grp)

  #browser()
  pps <- pds  %>%
    group_by(!!grps,bin, .draw) %>%
    summarise(tot = sum(tot_by_bin), .pred = sum(.prediction)) %>%
    group_by(!!grps, .draw) %>%
    mutate( .pred = (.pred/sum(.pred))*sum(tot)) %>%
    group_by(!!grps,bin, tot) %>%
    median_qi(.pred)

  p <- ggplot(pps) +
    geom_pointrange(aes(x=bin,y=.pred,ymin=.lower,ymax=.upper)) +
    geom_point(aes(x=bin,y=tot), col='skyblue') +
    xlab(xlab) +
    ylab('Prediction') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))

  if(!is.null(grp)) p <- p+facet_wrap(vars(!!grps),scales = 'free_y')

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
#' @param grp group to plot effects for, needs to include cvar if present
#' @param form Unused; retained for backwards compatibility.
#' @param grid optional plotting with facet grid along two variables (need to be in grp)
#' @param cvar continuous variable if used in model
#' @return A ggplot of median multipliers and intervals per bin, faceted by
#'   the first variable in \code{grp} (or by \code{grid}).
#'
#' @import cowplot
#' @import dplyr
#' @importFrom posterior draws_df
#' @export
#'
Fx_plot <- function(mod, grp='Year', form='~(1|bin:Year)', grid=NULL, cvar=NULL){

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

  if(class(mod$model)=='brmsfit'){

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

  } else if(class(mod$model)=='sdmTMB'){

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
    xlab('Length bin (cm)') +
    ylab('Multiplier') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))
  #browser()
  if(is.null(grid))  p <- p+facet_wrap(facets = vars(!!grps[[1]]),scales = 'free_y')
  if(!is.null(grid)) p <- p+facet_grid(rows = vars(!!as.symbol(grid[1])),cols = vars(!!as.symbol(grid[2])),scales = 'free_y')

  return(p)
}

#' Scale predicted compositions to catch
#'
#' Predicts the composition for each stratum in \code{scale_df} from a fitted
#' model and multiplies the predicted proportions by the stratum weight
#' (\code{predvar}, e.g. catch or predicted CPUE). Strata are the unique combinations of \code{grps} (plus
#' the coordinates for spatial models). The result is a set of draws of
#' numbers-at-bin per stratum, which can be summed over any subset of
#' \code{grps} and normalised to obtain weighted compositions (e.g. with
#' \code{scaled_comp_plot()}).
#'
#' The two settings of \code{remove_nuisance} correspond to two uses:
#' \describe{
#'   \item{Catch weighting (\code{remove_nuisance = FALSE}, default)}{The
#'     aim is the composition of the total catch of the full fleet. All model
#'     terms are evaluated at the values of each stratum in \code{scale_df},
#'     including the nuisance random effects (\code{re_form}). Every nuisance
#'     variable must therefore be in \code{grps}. \code{predvar} is the catch.}
#'   \item{CPUE-based scaling (\code{remove_nuisance = TRUE})}{The aim is a
#'     composition representative of the population, standardised for
#'     variables that affect the local availability of sizes to the sampling
#'     (e.g. fleet, season). Nuisance random effects are removed (population
#'     level); fixed effects, smooths and spatial fields are retained.
#'     \code{predvar} is typically predicted CPUE or abundance per stratum.}
#' }
#'
#' In both cases:
#' \itemize{
#'   \item Fixed effects (\code{ffx_form}) are evaluated at the levels in
#'     \code{scale_df}; these variables must therefore be in \code{grps}.
#'   \item Continuous covariates (\code{cvars}) are evaluated at their mean
#'     within each stratum of \code{scale_df}. Setting them to a reference
#'     value in \code{scale_df} standardises the compositions to that value.
#'   \item Spatial and spatio-temporal fields (spatial sdmTMB fits) are
#'     evaluated at the coordinates in \code{scale_df}.
#' }
#'
#' Nuisance levels that were not in the fitted data (e.g. unsampled vessels)
#' are handled differently by the backends. brms uses, in each draw, the
#' posterior draw of a randomly chosen fitted level
#' (\code{sample_new_levels = "uncertainty"}). sdmTMB predicts such rows at
#' the population level of all nuisance random effects, and a warning is
#' issued.
#'
#' @param scale_df data frame used for scaling, with one or more rows per
#'   stratum, containing \code{grps}, \code{predvar}, any \code{cvars}, and
#'   (spatial fits) the coordinate and time columns.
#' @param predvar Column in \code{scale_df} used for scaling (e.g. catch)
#' @param grps Grouping variables defining the scaling strata. Must include the
#'   \code{ffx_form} variables and, for catch weighting, the \code{re_form}
#'   variables.
#' @param fit A \code{ComPoM} object from \code{fit_model()}
#' @param iters Number of posterior draws (brms only; default uses all stored draws)
#' @param pgrid Reserved for future use
#' @param remove_nuisance Logical. \code{FALSE} (default) selects catch
#'   weighting, with nuisance random effects evaluated at each stratum's level.
#'   \code{TRUE} selects CPUE-based scaling, with nuisance random effects
#'   removed.
#' @return A data frame of draws (\code{.draw}) with predicted numbers-at-bin
#'   (\code{tot_by_bin}) per stratum and bin. The weight per stratum (sum of
#'   \code{predvar}) is in \code{n}.
#'
#' @import cowplot
#' @export
#'
scale_comps <- function(scale_df,
                        predvar    = 'catch',
                        fit        = NULL,
                        grps       = NULL,
                        iters      = NULL,
                        pgrid      = NULL,
                        remove_nuisance = FALSE) {

  if (is.null(fit)) stop("'fit' must be a ComPoM object from fit_model()")

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

  # ── nuisance random effects (re_form variables) ─────────────────────────────
  # Catch weighting (remove_nuisance = FALSE): the composition of the total catch
  # is predicted, so every nuisance variable is evaluated at the stratum's own
  # level. All of them must therefore be in grps (and scale_df).
  # CPUE-based scaling (remove_nuisance = TRUE): nuisance effects are removed, so
  # the composition is standardised to the population level of those variables.
  re_vars <- fit$re_vars
  if (!remove_nuisance && length(re_vars) > 0) {
    missing_re <- setdiff(re_vars, names(nd))
    if (length(missing_re) > 0)
      stop("Catch weighting (remove_nuisance = FALSE) predicts each stratum at its own ",
           "level of every nuisance variable; add ",
           paste0("'", missing_re, "'", collapse = ", "),
           " to `grps` and `scale_df`, or use remove_nuisance = TRUE.")
  }
  keep_re <- !remove_nuisance && length(re_vars) > 0

  if (inherits(fit$model, "brmsfit")) {

    # Levels not seen in the fit (e.g. unsampled vessels) take, in each draw,
    # the posterior draw of a randomly chosen fitted level
    # (brms sample_new_levels = "uncertainty").
    return(tidybayes::add_predicted_draws(
      nd, fit$model,
      allow_new_levels = TRUE,
      ndraws           = iters,
      value            = "tot_by_bin",
      re_formula       = if (keep_re) NULL else NA
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
#' @param scales passed to \code{facet_wrap()}/\code{facet_grid()}
#' @param comp_are,comp_are2 Optional data frames with \code{grps}, \code{bin}
#'   and \code{cvar}, overlaid as blue and orange points respectively.
#' @param cvar Column in \code{comp_are}/\code{comp_are2} holding the counts
#'   or catch per bin.
#' @return A ggplot.
#'
#' @import cowplot
#' @export
#'
scaled_comp_plot <- function(scaled_comp=NULL,
                             grps = NULL,
                             scaled=T,
                             scales = "free_y",
                             comp_are = NULL,
                             comp_are2 = NULL,
                             cvar='tcatch'){

  if(scaled == T) {

    preda <- scaled_comp %>%
      group_by(across(all_of(grps) ), bin, .draw) %>%
      summarise(tot = sum(tot_by_bin)) %>%
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
    {if(!is.null(comp_are)) geom_point(aes(x=bin,y=prop), col='skyblue',data=comp)} +
    {if(!is.null(comp_are2)) geom_point(aes(x=bin,y=prop), col='orange',data=comp2)} +
    xlab('Length bin (cm)') +
    ylab('Prediction') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=45, hjust=1,size = 8))

  if(length(grps)==1) p <- p+facet_wrap(facets = vars(!!as.symbol(grps)),scales = scales,drop = T,ncol = 4)
  if(length(grps)>=2) p <- p+facet_grid(rows = vars(!!as.symbol(grps[1])),cols = vars(!!!syms(grps[2:length(grps)])),scales = scales)

  p
}


#' Plot scaled compositions as ridgeline plot
#' @param scaled_comp data frame from \code{scale_comps()}
#' @param grps character vector of grouping variables used in \code{scale_comps()}
#' @param plot_y variable name (string) to stack on the y-axis (e.g. \code{"Year"})
#' @param plot_facet optional variable name (string) to facet by. \code{NULL} = no facet.
#' @param scales passed to \code{facet_wrap()}, default \code{"free_y"}
#' @param bin_range optional numeric vector \code{c(min, max)} to restrict plotted bins
#'
#' @import cowplot
#' @import ggridges
#' @export
#'
scaled_ridge_plot <- function(scaled_comp = NULL,
                              grps        = NULL,
                              plot_y      = NULL,
                              plot_facet  = NULL,
                              scales      = "free_y",
                              bin_range   = NULL) {

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
      fill      = "steelblue",
      alpha     = 0.6,
      color     = "grey30",
      linewidth = 0.3
    ) +
    ggplot2::xlab("Length bin (cm)") +
    ggplot2::ylab(plot_y) +
    ggplot2::theme_bw()

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
#' @export
#' @import ggplot2
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
                         fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(colour = "steelblue", linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::facet_wrap(~bin, labeller = ggplot2::label_both, scales = "free_y") +
    ggplot2::labs(
      x     = cvar,
      y     = "Multiplier (response scale, centred)",
      title = paste0("Marginal smooth effects by bin  -  ", cvar)
    ) +
    ggplot2::theme_bw(base_size = 10)
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
#' @export
#' @import ggplot2
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
    ggplot2::scale_fill_distiller(
      palette = "RdBu", limits = c(-lim, lim),
      name    = "Spatial\ndeviation\n(log scale)"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x     = lon_col,
      y     = lat_col,
      title = "Spatial field by bin (deviation from bin mean)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
