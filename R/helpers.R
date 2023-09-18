#' Selectivity
#' @param bins Bin structure
#' @param L50 50 pc selectivity
#' @param l95 95 pc selectivity offset from L50
#'
#' @export
#'
Selectivity <-function (bins,
                        L50,
                        L95) {

  Selex=1/(1+exp(-log(19)*( (bins-L50)/(L95-L50) )))

  return(Selex)

}
