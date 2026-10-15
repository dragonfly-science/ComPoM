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
"_PACKAGE"
