#' Title
#'
#' @param TSL data.table. Time Series List 
#' @param xi numeric. Damping ratio
#' @param Units character. Units of the time series
#'
#' @return
#' @export
#'
#' @examples
build_SDOF<- function(TSL,xi=0.05,Units){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  . <- NULL
  PSW <- TSL[ID=="AT",get_Spectra(.x=.SD,Units = Units,xi=xi),by=.(RecordSN,OCID,DIR)]
  
}
