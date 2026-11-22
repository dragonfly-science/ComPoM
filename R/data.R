#' Simulated example data
#'
#' Simulated length-composition sampling of a fishery, used in the examples and
#' the README. Fish get larger over the years, vessels select different sizes,
#' winter fish are larger, fish size varies along a spatial gradient, and sea-surface
#' temperature shifts the
#' composition. \code{ex_comp} holds the sampled length frequencies and
#' \code{ex_catch} the catch of the strata, for use with \code{scale_comps()}.
#'
#' @format \code{ex_comp} is a data frame with one row per stratum and length
#'   bin (144 strata x 6 bins):
#' \describe{
#'   \item{Year}{Factor, 2015 to 2020.}
#'   \item{Season}{Factor, \code{"summer"} or \code{"winter"}.}
#'   \item{Vessel}{Factor, \code{"V1"} to \code{"V4"}.}
#'   \item{trip}{Trip number within year, season and vessel.}
#'   \item{x, y}{Sampling location (arbitrary projected coordinates, 0 to 100).}
#'   \item{sst}{Sea-surface temperature at the sample.}
#'   \item{length}{Length bin (10 to 60).}
#'   \item{count}{Number of fish measured in the bin (60 fish per stratum).}
#' }
#'
#' \code{ex_catch} is a data frame with one row per stratum (144 rows) and the
#' stratum variables of \code{ex_comp} (\code{Year}, \code{Season},
#' \code{Vessel}, \code{trip}, \code{x}, \code{y}, \code{sst}) plus
#' \code{catch}, the total catch of the stratum.
#' @source Simulated, see \code{data-raw/ex_data.R}.
#' @name ex_data
#' @aliases ex_comp ex_catch
#' @examples
#' head(ex_comp)
#' head(ex_catch)
NULL
