#' Title
#'
#' @param .x data.table
#' @param s vector
#' @param t vector
#' @param dt double
#' @param Units string
#' @param TargetUnits string
#' @param Tn vector
#' @param xi double
#'
#' @return list
#' @export 
#'
#' @import expm
#' @import data.table
#' @importFrom MASS pinv
#' @examples
get_Spectra <- function(.x=NULL,s=NULL,t=NULL,dt=NULL,Tn = NULL,xi=0.05){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  # Check internal  -------------------------------------------------------------
  if(is.null(s) && is.null(.x)){
    stop("No acceleration data provided.")
  }
  if(!is.null(.x) & is.data.table(.x)){
    AT <- .x$s
    t <- .x$t
  }
  if(is.null(t) && is.null(.x) && is.null(dt) ){
    stop("No time data provided.")
  }
  if(is.null(.x) && !is.null(s)){
    AT <-  copy(s)
  }
  if(!is.null(t)){
    dt <- diff(t) |> mean()
  }
  if(is.null(Tn)){
    seq1 <- (10^seq(log10(0.001), log10(0.1), length.out = 40)) |> round(3)
    seq2 <- (10^seq(log10(0.1), log10(1), length.out = 20)) |> round(2)
    seq3 <- (10^seq(log10(1), log10(10), length.out = 10)) |> round(1)
    Tn <- c(0,seq1, seq2, seq3) |> unique() 
  }
  
  PGA <- max(abs(AT))
  PSA <- numeric(length(Tn))
  PSA <- sapply(Tn, function(Ti) {
    if(Ti>0){
      Wn <- 2 * pi / Ti
      A <- matrix(c(0, 1, -Wn^2, -2 * xi * Wn), nrow = 2, byrow = TRUE)
      Ae <- expm(A * dt)
      AeB <- (Ae - diag(2)) %*% MASS::ginv(A) %*% matrix(c(0, 1), nrow = 2)
      
      y <- matrix(0, nrow = 2, ncol = length(AT))
      for (k in 2:length(AT)) {
        y[, k] <- Ae %*% y[, k-1] + AeB * AT[k]
      }
      return(max(abs(y[1, ])) * Wn^2)
    } else {
      return(PGA)
    }
    
  })
  ## Pack Spectra  ----

  return(list( Tn=Tn, xi=xi,PSA = PSA))

}
