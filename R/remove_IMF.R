
#' Title
#'
#' @param X vector. Time Series
#' @param EMD vector. Empirical Mode Decomposition of X
#' @param removeIMF1 integer. Number of IMFs to remove from the beginning
#' @param removeIMFn integer. Number of IMFs to remove from the end
#'
#' @return list
#' @export
#'
#' @examples
remove_IMF <- function(X,EMD,removeIMF1=0,removeIMFn=0){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  if (removeIMF1>0 ||removeIMFn>0) {
    nimf <- EMD$nimf
    i <- removeIMF1
    j <- removeIMFn
    COLS <- colnames(EMD$imf)[(i+1):(nimf-j)]
    X <-  EMD$imf[,COLS,with = FALSE] |> rowSums()
  }
  return(X)
}
