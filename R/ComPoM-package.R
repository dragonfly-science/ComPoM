#' ComPoM: Poisson-multinomial models to standardise and scale compositions
#'
#' ComPoM fits Poisson-multinomial models to compositional count data (e.g.
#' length or age frequencies) and uses them to predict compositions for strata
#' with known catch. The predictions are aggregated into catch-weighted
#' compositions with uncertainty.
#'
#' The workflow consists of four steps:
#' \enumerate{
#'   \item \code{data_prep()}: aggregate counts to strata x bins, fill zero
#'     counts and add stratum totals.
#'   \item \code{fit_model()}: fit the model with \code{brms} or
#'     \code{sdmTMB}, optionally with smooths of continuous covariates and
#'     bin-specific spatial and spatio-temporal random fields.
#'   \item Diagnostics and effects: \code{post_pred_group()},
#'     \code{Fx_plot()}, \code{cvar_smooth_plot()},
#'     \code{spatial_field_plot()}.
#'   \item Scaling: \code{scale_comps()}, then \code{scaled_comp_plot()} or
#'     \code{scaled_ridge_plot()}.
#' }
#'
#' See \code{vignette("ComPoM")} for the model equations and an example.
#'
#' @keywords internal
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import rlang
#' @import cowplot
#' @importFrom tidybayes mean_qi
#' @importFrom stats quantile setNames
"_PACKAGE"

# Column names used in tidy evaluation (dplyr, ggplot2) and by the posterior draw helpers
utils::globalVariables(c(
  ".", ".draw", ".linpred", ".lower", ".pad", ".pred", ".prediction", ".upper",
  "bin", "est", "est_c", "hi", "iter", "lo", "lp", "ml", "pred", "prop", "tot",
  "tot_by_bin", "var", "level", "sign", "share", "effect", "lo", "hi", "p", "p_r", "S", "S_r",
  "dS", "infl", "w", "Lower", "Upper", "panel", "role", "component", "step", "x", "m", "v",
  "j", "faded", "b"
))
