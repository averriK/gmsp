#' Title
#'
#' @param ATL data.table
#' @param TargetUnits character
#' @param kh vector
#' @param full boolean
#'
#' @return
#' @export
#'
#' @examples



get_Newmark <- function(ATL,TargetUnits="mm",kh=c(0.01,0.02,0.05, 0.10, 0.15, 0.20, 0.25, 0.30,0.35, 0.40,0.45, 0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50,3.00),full=FALSE){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  . <- NULL
  g <- .getG(TargetUnits) #GMSP$g
  # Newmark Displacements  -------------------------------------------------------------
  
  NDL <- rbindlist(
    sapply(kh, function(x) {
      ATL[, .(
        ID = sprintf("DN%02d", round(x * 100)), 
        value = get_ND(AT = s, t = t, kh = x, FULL = full)), 
        by = .(RecordSN, DIR, OCID)]}, simplify = FALSE))
  NDW <- dcast(NDL, RecordSN + OCID + DIR ~ ID, value.var = "value")
  return(list(NDL=NDL,NDW=NDW))
}

