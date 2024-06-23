#' Title
#'
#' @param TSL data.table. Time Series List 
#' @param xi numeric. Damping ratio
#'
#' @return
#' @export
#'
#' @examples
build_SDOF<- function(TSL,xi=0.05){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  PSW <- TSL[ID=="AT",get_Spectra(.x=.SD,xi=xi),by=list(RecordSN,OCID,DIR)]
  
}
