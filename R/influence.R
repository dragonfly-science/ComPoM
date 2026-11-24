#' Statistics of a composition
#'
#' Helpers for the \code{stat} argument of \code{step_comps()} and
#' \code{influ_comps()}. A statistic is a function of the bin values and the
#' proportions in each bin that returns one number.
#'
#' @param prob Probability of the quantile, between 0 and 1.
#' @return A function with arguments \code{bin} and \code{p}.
#' @examples
#' comp_quantile(0.9)(bin = c(10, 20, 30), p = c(0.2, 0.5, 0.3))
#' @export
comp_quantile <- function(prob) {
  force(prob)
  function(bin, p) {
    o  <- order(bin)
    cp <- cumsum(p[o]) / sum(p)
    bin[o][which(cp >= prob)[1]]
  }
}

# "mean", a probability in (0, 1) for a quantile, or a function(bin, p)
comp_stat_fun <- function(stat) {
  if (is.function(stat)) return(stat)
  if (identical(stat, "mean")) return(function(bin, p) sum(bin * p) / sum(p))
  if (is.numeric(stat) && length(stat) == 1 && stat > 0 && stat < 1) return(comp_quantile(stat))
  stop("`stat` must be \"mean\", a probability between 0 and 1 (quantile), or a function(bin, p)")
}

comp_stat_label <- function(stat) {
  if (is.function(stat)) "Statistic"
  else if (identical(stat, "mean")) "Mean bin"
  else paste0(round(100 * stat), "th percentile of bin")
}

# Draws of numbers-at-bin summed to `by`; statistic per draw and level of `by`
stat_draws <- function(sc, by, stat) {
  f <- comp_stat_fun(stat)
  sc |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "bin", ".draw")))) |>
    dplyr::summarise(v = sum(tot_by_bin), .groups = "drop") |>
    dplyr::mutate(bin = as.numeric(as.character(bin))) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, ".draw")))) |>
    dplyr::summarise(S = f(bin, v), .groups = "drop")
}

# Draws of proportions per level of `by` and bin
prop_draws <- function(sc, by) {
  sc |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "bin", ".draw")))) |>
    dplyr::summarise(v = sum(tot_by_bin), .groups = "drop") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(by, ".draw")))) |>
    dplyr::mutate(p = v / sum(v)) |>
    dplyr::ungroup() |>
    dplyr::mutate(bin = as.numeric(as.character(bin)))
}

summ_draws <- function(df, value, by) {
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::summarise(mean  = mean(.data[[value]]),
                     Lower = stats::quantile(.data[[value]], 0.025),
                     Upper = stats::quantile(.data[[value]], 0.975),
                     .groups = "drop")
}

# Predictions of one component, with an optional set of extra terms set to a reference
predict_component <- function(fit, component, scale_df, predvar, grps, grid, selectivity,
                              reference, weight, density, extra = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  terms_in_model <- unique(c(fit$ffx_vars, fit$cvars, fit$re_vars))
  cat_all <- intersect(union(selectivity, extra), terms_in_model)
  if (component == "catch") {
    suppressWarnings(scale_comps(scale_df, predvar = predvar, fit = fit, grps = grps,
                                 mode = "catch", selectivity = extra, reference = reference))
  } else {
    suppressWarnings(scale_comps(fit = fit, mode = "index", grps = intersect(grps, names(grid)),
                                 grid = grid,
                                 selectivity = cat_all, reference = reference,
                                 weight = weight, density = intersect(density, terms_in_model)))
  }
}

#' Add model terms one at a time and track a statistic of the composition
#'
#' Refits the model with terms added successively (a step plot, as used for
#' CPUE standardisation) and computes, for each step, a statistic of the
#' annual composition (by default the mean bin, e.g. mean length) with
#' uncertainty. The statistic is computed for the catch composition
#' (\code{mode = "catch"} in \code{scale_comps()}) and for the index
#' composition (\code{mode = "index"}).
#'
#' The first step is the model with bin-specific intercepts only. Terms are
#' added in the order of the model formula (fixed-effect variables, then
#' random-effect variables, then covariates, then the spatial fields), or in
#' the order given by \code{order}. Only the \code{"TMB"} backend is supported.
#'
#' @param fit A \code{ComPoM} object from \code{fit_model()} with
#'   \code{backend = "TMB"}.
#' @param scale_df,predvar Data frame and weight column for the catch
#'   component; see \code{scale_comps()}.
#' @param grps Grouping variables for the scaling; see \code{scale_comps()}.
#' @param by Variable in \code{grps} (or \code{grid}) over which the statistic
#'   is plotted, usually the year.
#' @param grid,selectivity,reference,weight,density Settings of the index
#'   component; see \code{scale_comps()}.
#' @param order Optional character vector with the terms in the order in which
#'   they are added. Terms are the \code{ffx_form} and \code{re_form}
#'   variables, the covariates, and \code{"spatial"} for spatial fits.
#' @param stat \code{"mean"} (default), a probability between 0 and 1 for a
#'   quantile of the bin (e.g. \code{0.9}), or a function of \code{bin} and
#'   \code{p}. See \code{comp_quantile()}.
#' @param components Components to compute: \code{"catch"} and/or
#'   \code{"index"}. A component is skipped when its data frame is not
#'   supplied.
#' @return A list of class \code{comp_step} with \code{steps} (statistic by
#'   step, component and level of \code{by}, with a 95\% interval),
#'   \code{summary} (degrees of freedom, log-likelihood and AIC per step),
#'   \code{by} and \code{stat_label}.
#' @seealso \code{plot_step_comps()}, \code{influ_comps()}
#' @examples
#' \donttest{
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year", "Season + Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2)
#' st <- step_comps(fit, scale_df = ex_catch, predvar = "catch",
#'                  grps = c("Year", "Season", "Vessel"), by = "Year",
#'                  grid = data.frame(Year = factor(2015:2020)))
#' plot_step_comps(st)
#' }
#' @export
step_comps <- function(fit, scale_df = NULL, predvar = "catch", grps, by,
                       grid = NULL, selectivity = NULL, reference = NULL,
                       weight = NULL, density = NULL, order = NULL, stat = "mean",
                       components = c("catch", "index")) {
  sp <- fit$spec
  if (is.null(sp)) stop("`fit` has no stored specification; refit with the current fit_model()")
  if (sp$backend != "TMB") stop("step_comps() supports the TMB backend only")
  comp_stat_fun(stat)
  components <- match.arg(components, c("catch", "index"), several.ok = TRUE)
  if ("catch" %in% components && is.null(scale_df)) components <- setdiff(components, "catch")
  if ("index" %in% components && is.null(grid)) components <- setdiff(components, "index")
  if (length(components) == 0) stop("Supply `scale_df` (catch) and/or `grid` (index)")

  ffx <- fit$ffx_vars
  res <- fit$re_vars
  cv  <- fit$cvars
  all_terms <- c(ffx, res, cv, if (!is.null(sp$mesh)) "spatial")
  if (is.null(order)) order <- all_terms
  if (!setequal(order, all_terms) || anyDuplicated(order))
    stop("`order` must contain each of: ", paste(all_terms, collapse = ", "))

  nz <- function(x) if (length(x) == 0) NULL else paste(x, collapse = " + ")
  out <- list(); summ <- list()
  for (k in 0:length(order)) {
    on  <- order[seq_len(k)]
    lab <- if (k == 0) "bin only" else paste("+", order[k])
    spatial_on <- "spatial" %in% on
    fit_k <- suppressMessages(suppressWarnings(fit_model(
      ffx_form = nz(intersect(ffx, on)), re_form = nz(intersect(res, on)),
      cvars = if (length(intersect(cv, on)) > 0) intersect(cv, on) else NULL,
      spl = sp$spl, knots = sp$knots, data = sp$data, backend = "TMB",
      mesh = if (spatial_on) sp$mesh else NULL,
      time = if (spatial_on) sp$time else NULL, coords = if (spatial_on) sp$coords else NULL,
      spatiotemporal = sp$spatiotemporal, share_spatial_sd = sp$share_spatial_sd,
      share_spatiotemporal_sd = sp$share_spatiotemporal_sd, priors = sp$priors,
      newton_loops = sp$newton_loops, iter = sp$iter, warmup = sp$warmup, thin = sp$thin,
      dist = sp$dist %||% "poisson"
    )))
    for (cmp in components) {
      sc <- predict_component(fit_k, cmp, scale_df, predvar, grps, grid, selectivity,
                              reference, weight, density)
      s <- summ_draws(stat_draws(sc, by, stat), "S", by)
      names(s)[names(s) == by] <- "level"
      s$component <- cmp; s$step <- lab
      out[[length(out) + 1]] <- s
    }
    summ[[k + 1]] <- data.frame(
      step = lab, df = length(fit_k$model$sd_report$par.fixed),
      logLik = as.numeric(stats::logLik(fit_k$model)), AIC = stats::AIC(fit_k$model),
      pdHess = isTRUE(fit_k$model$sd_report$pdHess))
  }
  steps <- dplyr::bind_rows(out)
  lv <- c("bin only", paste("+", order))
  steps$step <- factor(steps$step, levels = lv)
  summary <- dplyr::bind_rows(summ)
  summary$step <- factor(summary$step, levels = lv)
  structure(list(steps = steps, summary = summary, by = by, stat_label = comp_stat_label(stat)),
            class = "comp_step")
}

#' Plot the step results
#'
#' One panel per step, with the statistic by level of \code{by} for the current
#' step (orange, with its 95\% interval), the previous step (skyblue) and
#' earlier steps (grey). Rows are the components (catch and index).
#'
#' @param step Output of \code{step_comps()}.
#' @param xlab,ylab Axis labels. \code{ylab} defaults to the statistic.
#' @return A ggplot.
#' @seealso \code{step_comps()}
#' @export
plot_step_comps <- function(step, xlab = step$by, ylab = step$stat_label) {
  d <- step$steps
  lv <- levels(d$step)
  lev_num <- suppressWarnings(as.numeric(as.character(d$level)))
  d$x <- if (!anyNA(lev_num)) lev_num else as.numeric(factor(d$level))
  bg <- dplyr::bind_rows(lapply(seq_along(lv), function(i) {
    d |>
      dplyr::filter(as.integer(step) <= i) |>
      dplyr::mutate(panel = factor(lv[i], levels = lv),
                    role = dplyr::case_when(as.integer(step) == i ~ "Current",
                                            as.integer(step) == i - 1 ~ "Previous",
                                            TRUE ~ "Earlier"))
  }))
  bg$role <- factor(bg$role, levels = c("Earlier", "Previous", "Current"))
  cur <- dplyr::filter(bg, role == "Current")
  ggplot(bg, aes(x = x, y = mean, group = step, colour = role)) +
    geom_ribbon(data = cur, aes(ymin = Lower, ymax = Upper), fill = "orange",
                colour = NA, alpha = 0.25) +
    geom_line(aes(linewidth = role)) +
    geom_point(data = cur, size = 1.2) +
    scale_colour_manual(values = c(Earlier = "grey65", Previous = "skyblue", Current = "orange"),
                        name = NULL) +
    scale_linewidth_manual(values = c(Earlier = 0.3, Previous = 0.6, Current = 0.8), guide = "none") +
    facet_grid(component ~ panel, scales = "free_y") +
    labs(x = xlab, y = ylab) +
    theme_compom(font_size = 10) +
    theme(legend.position = "bottom", strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Influence of model terms on the composition
#'
#' For each term, computes (1) the effect of each level on every bin, (2) the
#' distribution of the data over the levels by year, and (3) the influence of
#' the term on the annual composition. The influence is the log ratio of the
#' annual proportion in each bin predicted with the term at the observed levels
#' to the proportion predicted with the term at its reference. The same
#' quantity is also given for a statistic of the composition (by default the
#' mean bin). The design follows the coefficient-distribution-influence plots
#' used for CPUE standardisation.
#'
#' The influence is computed for the catch component
#' (\code{mode = "catch"} in \code{scale_comps()}) and the index component
#' (\code{mode = "index"}). Terms in \code{selectivity} have no influence on the
#' index by construction and are given zero influence for that component.
#' Random effects are removed with the TMB backend only together, so that with
#' more than one random-effect variable they form one term (labelled with all
#' variables). Effects of covariate smooths and spatial fields are not
#' included; see \code{cvar_smooth_plot()} and \code{spatial_field_plot()}.
#'
#' @inheritParams step_comps
#' @param terms Terms to evaluate. Default: the \code{ffx_form} and
#'   \code{re_form} variables except \code{by}.
#' @param components Components to compute: \code{"catch"} and/or
#'   \code{"index"}.
#' @return A list of class \code{comp_influ} with data frames \code{effects}
#'   (log multiplier by term, level and bin), \code{distribution} (share of the
#'   weight by term, level, component and year), \code{influence} (log ratio by
#'   term, component, year and bin) and \code{influence_stat} (difference in
#'   the statistic by term, component and year).
#' @seealso \code{plot_cdi_comps()}, \code{step_comps()}
#' @examples
#' \donttest{
#' dat <- data_prep(ex_comp, c("Year", "Season", "Vessel", "trip"), "length", "count")
#' fit <- fit_model("Year + Season", "Vessel", data = dat, backend = "TMB",
#'                  iter = 300, warmup = 100, thin = 2)
#' inf <- influ_comps(fit, scale_df = ex_catch, predvar = "catch",
#'                    grps = c("Year", "Season", "Vessel"), by = "Year",
#'                    terms = "Season")
#' plot_cdi_comps(inf, "Season")
#' }
#' @export
influ_comps <- function(fit, scale_df = NULL, predvar = "catch", grps, by,
                        terms = NULL, grid = NULL, selectivity = NULL, reference = NULL,
                        weight = NULL, density = NULL, stat = "mean",
                        components = c("catch", "index")) {
  if (!inherits(fit$model, "sdmTMB")) stop("influ_comps() supports the TMB backend only")
  comp_stat_fun(stat)
  components <- match.arg(components, c("catch", "index"), several.ok = TRUE)
  if ("catch" %in% components && is.null(scale_df)) components <- setdiff(components, "catch")
  if ("index" %in% components && is.null(grid)) components <- setdiff(components, "index")
  if (length(components) == 0) stop("Supply `scale_df` (catch) and/or `grid` (index)")

  if (is.null(terms)) terms <- setdiff(c(fit$ffx_vars, fit$re_vars), by)
  if (by %in% terms) stop("`by` (", by, ") is the axis of the influence and cannot be a term")
  bad <- setdiff(terms, c(fit$ffx_vars, fit$re_vars))
  if (length(bad) > 0) stop("Terms not supported (fixed-effect or random-effect variables only): ",
                            paste(bad, collapse = ", "))
  # with TMB the random effects can only be removed together
  re_t <- intersect(terms, fit$re_vars)
  groups <- as.list(setdiff(terms, re_t)); names(groups) <- unlist(groups)
  if (length(re_t) > 0) {
    lab <- paste(fit$re_vars, collapse = " + ")
    if (length(fit$re_vars) == 1) lab <- fit$re_vars
    groups[[lab]] <- fit$re_vars
  }

  seed <- sample.int(1e6, 1)
  eff <- list(); dist <- list(); infl <- list(); infs <- list()
  for (cmp in components) {
    base <- predict_component(fit, cmp, scale_df, predvar, grps, grid, selectivity,
                              reference, weight, density, seed = seed)
    pb <- prop_draws(base, by); sb <- stat_draws(base, by, stat)
    wdf <- if (cmp == "catch") scale_df else {
      g <- as.data.frame(grid); g$.w <- if (is.null(weight)) 1 else g[[weight]]; g
    }
    wcol <- if (cmp == "catch") predvar else ".w"
    for (lab in names(groups)) {
      vars <- groups[[lab]]
      if (cmp == "index" && all(vars %in% c(selectivity, setdiff(fit$re_vars, density)))) {
        # selectivity terms have no influence on the index by construction
        lv <- unique(pb[[by]])
        z <- expand.grid(level = lv, bin = unique(pb$bin))
        z$mean <- 0; z$Lower <- 0; z$Upper <- 0
        names(z)[1] <- by
        z$term <- lab; z$component <- cmp
        infl[[length(infl) + 1]] <- z
        zs <- data.frame(level = lv, mean = 0, Lower = 0, Upper = 0); names(zs)[1] <- by
        zs$term <- lab; zs$component <- cmp
        infs[[length(infs) + 1]] <- zs
      } else {
        red <- predict_component(fit, cmp, scale_df, predvar, grps, grid, selectivity,
                                 reference, weight, density, extra = vars, seed = seed)
        pr <- prop_draws(red, by); sr <- stat_draws(red, by, stat)
        j <- dplyr::inner_join(pb, pr, by = c(by, "bin", ".draw"), suffix = c("", "_r")) |>
          dplyr::mutate(infl = log(p / p_r))
        s <- summ_draws(j, "infl", c(by, "bin")); s$term <- lab; s$component <- cmp
        infl[[length(infl) + 1]] <- s
        js <- dplyr::inner_join(sb, sr, by = c(by, ".draw"), suffix = c("", "_r")) |>
          dplyr::mutate(dS = S - S_r)
        ss <- summ_draws(js, "dS", by); ss$term <- lab; ss$component <- cmp
        infs[[length(infs) + 1]] <- ss
      }
      # distribution of the weight over the levels of each variable by year
      for (v in vars) {
        if (v %in% names(wdf) && by %in% names(wdf)) {
          d <- wdf |>
            dplyr::group_by(dplyr::across(dplyr::all_of(unique(c(by, v))))) |>
            dplyr::summarise(w = sum(.data[[wcol]]), .groups = "drop") |>
            dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
            dplyr::mutate(share = w / sum(w)) |>
            dplyr::ungroup()
          names(d)[names(d) == v] <- "level"
          d$level <- as.character(d$level); d$var <- v; d$term <- lab; d$component <- cmp
          dist[[length(dist) + 1]] <- d
        }
      }
    }
  }
  # effect of each level on each bin (multiplier relative to the geometric mean over levels)
  for (v in unlist(groups)) {
    e <- fx_effects(fit, grp = v)
    e <- data.frame(var = v, level = as.character(e[[v]]), bin = as.numeric(as.character(e$bin)),
                    effect = log(e$pred), lo = log(e$.lower), hi = log(e$.upper))
    eff[[length(eff) + 1]] <- e
  }
  structure(list(effects = dplyr::bind_rows(eff), distribution = dplyr::bind_rows(dist),
                 influence = dplyr::bind_rows(infl), influence_stat = dplyr::bind_rows(infs),
                 groups = groups, by = by, stat_label = comp_stat_label(stat)),
            class = "comp_influ")
}

#' Coefficient, distribution and influence plot for one term
#'
#' A coefficient-distribution-influence plot (Bentley et al. 2011) for one term
#' of \code{influ_comps()}, with the axes shared between the panels.
#' \itemize{
#'   \item Centre: the distribution of the data. Bubbles show the share of the
#'     data in each level of the term (x) in each year (y).
#'   \item Top, sharing the level axis (x): the effect of each level on every bin
#'     (y), as the log of the multiplier relative to the geometric mean over
#'     levels. Bars extend to the right (more fish in the bin than average) or to
#'     the left (fewer), with 95\% intervals. Bar lengths are scaled to the
#'     largest effect, and a scale bar gives the reference length. With many
#'     levels, a heat map replaces the bars. Cells whose interval includes zero
#'     are faded.
#'   \item Right, sharing the year axis (y): the influence of the term on the
#'     annual composition in every bin (x), as the log ratio of the predicted
#'     proportion with and without the term. Cells whose interval includes zero
#'     are faded.
#' }
#' Orange is below and skyblue above the reference.
#'
#' @param influ Output of \code{influ_comps()}.
#' @param term Term (a name in \code{influ$groups}, i.e. a variable, or the
#'   joint random-effect label).
#' @param component \code{"catch"} or \code{"index"}, for the distribution and
#'   influence panels.
#' @param var Variable whose levels are shown, for a term that groups several
#'   variables. Default: the first.
#' @param style \code{"auto"} (bars for up to \code{max_levels} levels, otherwise a
#'   heat map), \code{"bars"} or \code{"tile"}.
#' @param max_levels Maximum number of levels shown as bars with
#'   \code{style = "auto"}.
#' @param xlab Label of the bin axis.
#' @return A patchwork of ggplots.
#' @references Bentley, N., Kendrick, T.H., Starr, P.J., Breen, P.A. (2012).
#'   Influence plots and metrics: tools for better understanding fisheries
#'   catch-per-unit-effort standardizations. ICES Journal of Marine Science 69,
#'   84-88.
#' @seealso \code{influ_comps()}
#' @export
plot_cdi_comps <- function(influ, term, component = NULL, var = NULL,
                           style = c("auto", "bars", "tile"), max_levels = 10,
                           xlab = "Length bin (cm)") {
  style <- match.arg(style)
  if (!term %in% names(influ$groups))
    stop("`term` must be one of: ", paste(names(influ$groups), collapse = ", "))
  vars <- influ$groups[[term]]
  v    <- if (is.null(var)) vars[1] else var
  if (!v %in% vars) stop("`var` must be one of: ", paste(vars, collapse = ", "))
  comp <- if (is.null(component)) unique(influ$influence$component)[1] else component
  by   <- influ$by
  neg <- "orange"; pos <- "skyblue"

  # levels (x) and years (y) as integer positions, so that panels share exactly the same axes
  e <- dplyr::filter(influ$effects, var == v)
  lvl_num <- suppressWarnings(as.numeric(e$level))
  lv <- unique(e$level)
  if (length(lv) > max_levels) {
    lv <- (e |> dplyr::group_by(level) |> dplyr::summarise(m = mean(effect)) |>
             dplyr::arrange(m))$level
  } else if (!anyNA(lvl_num)) {
    lv <- unique(e$level[order(lvl_num)])
  }
  nl <- length(lv)
  e$j <- match(e$level, lv)

  yrs <- unique(as.character(influ$influence[[by]]))
  yn  <- suppressWarnings(as.numeric(yrs))
  if (!anyNA(yn)) yrs <- yrs[order(yn)]
  ny  <- length(yrs)

  bins <- sort(unique(e$bin))
  bstep <- if (length(bins) > 1) min(diff(bins)) else 1
  xl <- c(0.4, nl + 0.6); yl <- c(0.4, ny + 0.6)
  xscale <- scale_x_continuous(breaks = seq_len(nl), labels = lv, limits = xl, expand = c(0, 0))
  yscale <- scale_y_continuous(breaks = seq_len(ny), labels = yrs, limits = yl, expand = c(0, 0))
  lab_theme <- if (nl > 12) theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7)) else theme()

  # centre: distribution of the data
  d <- dplyr::filter(influ$distribution, term == !!term, var == v, component == comp)
  d$j <- match(d$level, lv); d$i <- match(as.character(d[[by]]), yrs)
  pb <- ggplot(d, aes(x = j, y = i, size = share)) +
    geom_point(colour = "skyblue", alpha = 0.85) +
    scale_size_area(max_size = 7, name = "Share of\ndata") +
    xscale + yscale +
    labs(x = v, y = by) +
    theme_compom(font_size = 10) + lab_theme

  # top: effect of each level on each bin
  use_bars <- style == "bars" || (style == "auto" && nl <= max_levels)
  e$faded <- ifelse(e$lo <= 0 & e$hi >= 0, 0.45, 1)
  lim_e <- max(abs(e$effect))
  if (use_bars) {
    k <- 0.42 / lim_e
    e$sign <- ifelse(e$effect >= 0, pos, neg)
    pa <- ggplot(e) +
      geom_segment(aes(x = j, xend = j, y = min(bin) - bstep / 2, yend = max(bin) + bstep / 2),
                   colour = "grey75", linewidth = 0.3) +
      geom_rect(aes(xmin = j, xmax = j + k * effect, ymin = bin - 0.4 * bstep,
                    ymax = bin + 0.4 * bstep, fill = sign, alpha = faded)) +
      geom_segment(aes(x = pmax(j - 0.45, j + k * lo), xend = pmin(j + 0.45, j + k * hi), y = bin,
                       yend = bin), colour = "grey40", linewidth = 0.25, alpha = 0.7) +
      # scale bar: bar widths are scaled to the largest effect, so a reference length is shown
      geom_segment(data = data.frame(b = signif(lim_e / 2, 1)),
                   aes(x = 0.5, xend = 0.5 + k * b, y = min(bins) - 0.9 * bstep,
                       yend = min(bins) - 0.9 * bstep), linewidth = 1, colour = "grey30") +
      geom_text(data = data.frame(b = signif(lim_e / 2, 1)),
                aes(x = 0.5 + k * b / 2, y = min(bins) - 1.5 * bstep, label = paste0(b, " (log)")),
                size = 2.6, colour = "grey30") +
      scale_y_continuous(breaks = pretty(range(bins)),
                         expand = expansion(add = c(1.9 * bstep, 0.6 * bstep))) +
      scale_fill_identity() + scale_alpha_identity()
  } else {
    pa <- ggplot(e) +
      geom_tile(aes(x = j, y = bin, fill = effect, alpha = faded), height = bstep) +
      scale_fill_gradient2(low = neg, mid = "white", high = pos, midpoint = 0,
                           limits = c(-lim_e, lim_e), name = "Effect\n(log)") +
      scale_alpha_identity()
  }
  pa <- pa + xscale + labs(x = NULL, y = xlab) +
    theme_compom(font_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # right: influence of the term by year and bin
  i <- dplyr::filter(influ$influence, term == !!term, component == comp)
  i$i <- match(as.character(i[[by]]), yrs)
  i$faded <- ifelse(i$Lower <= 0 & i$Upper >= 0, 0.45, 1)
  lim_i <- max(abs(i$mean), 1e-8)
  pc <- ggplot(i, aes(x = bin, y = i, fill = mean, alpha = faded)) +
    geom_tile(width = bstep) +
    scale_fill_gradient2(low = neg, mid = "white", high = pos, midpoint = 0,
                         limits = c(-lim_i, lim_i), name = "Influence\n(log ratio)") +
    scale_alpha_identity() +
    yscale + labs(x = xlab, y = NULL) +
    theme_compom(font_size = 10) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin(t = 14, r = 5.5, b = 5.5, l = 5.5))

  design <- "AA#\nBBC"
  patchwork::wrap_plots(A = pa, B = pb, C = pc, design = design) +
    patchwork::plot_layout(widths = c(4, 4, 2.2), heights = c(2.2, 3.5)) +
    patchwork::plot_annotation(tag_levels = "a")
}
