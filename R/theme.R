#' ComPoM plot theme
#'
#' The \code{cowplot} theme with faded major gridlines, used by all plotting
#' functions in the package. Use it (or \code{ggplot2::theme_set(theme_compom())})
#' to match other figures to the package plots.
#'
#' @param font_size Base font size.
#' @return A ggplot2 theme.
#' @examples
#' ggplot2::ggplot(ex_comp, ggplot2::aes(length, count)) +
#'   ggplot2::geom_col() +
#'   theme_compom()
#' @export
theme_compom <- function(font_size = 11) {
  cowplot::theme_cowplot(font_size = font_size) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3))
}
