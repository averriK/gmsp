#' Title
#'
#' @param .x data.table
#' @param a vector
#' @param t vector
#' @param dt double
#' @param Units string
#' @param TargetUnits string
#' @param Tn vector
#' @param xi double
#'
#' @return list
#' @export buildPSA
#'
#' @import expm
#' @import data.table
#' @importFrom MASS pinv
#' @examples
buildPSA <- function(.x=NULL,a=NULL,t=NULL,dt=NULL,Units,TargetUnits="mm",
    Tn = c(0.01, 0.02, 0.022, 0.025, 0.029, 0.03, 0.032, 0.035, 0.036, 0.04, 0.042, 0.044, 0.045, 0.046, 0.048, 0.05, 0.055, 0.06, 0.065, 0.067, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.11, 0.12, 0.125, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.3, 0.32, 0.34, 0.35, 0.36, 0.38, 0.4, 0.42, 0.44, 0.45, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 0.667, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.5, 2.6, 2.8, 3, 3.2, 3.4, 3.5, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 7.5, 10),
    xi=0.05){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  # Check internal  -------------------------------------------------------------
  if(is.null(a) && is.null(.x)){
    stop("No acceleration data provided.")
  }
  if(!is.null(.x) & is.data.table(.x)){
    AT <- .x$s
    t <- .x$t
  }
  if(is.null(t) && is.null(.x) && is.null(dt) ){
    stop("No time data provided.")
  }
  if(is.null(.x) && !is.null(a)){
    AT <-  copy(a)
  }
  if(!is.null(t)){
    dt <- diff(t) |> mean()
  }
  if(grepl(Units,pattern = "[///+]")){
    Units <- (str_split(Units, pattern = "[///+]") |> unlist())[1]
  }
  if(!(tolower(Units) %in% c("mm","cm","m","gal","g"))){
   stop("Invalid Units")
  }
  if (tolower(Units) != TargetUnits) {
    SFU <- .getSF(SourceUnits = tolower(Units), TargetUnits = TargetUnits)
    # AT <- map(AT,function(x){x*SFU})
    # ATo[, (colnames(ATo)) := lapply(.SD, function(x) {x * SFU})]
    AT <- AT * SFU
  }

  PSA <- numeric(length(Tn))
  PSV <- numeric(length(Tn))
  SD <- numeric(length(Tn))


  for (n in seq_along(Tn)) {
    Wn <- 2 * pi / Tn[n]
    A <- matrix(c(0, 1, -Wn^2, -2 * xi * Wn), nrow = 2, byrow = TRUE)
    Ae <- expm::expm(A * dt)
    AeB <- (Ae - diag(2)) %*% MASS::ginv(A) %*% matrix(c(0, 1), nrow = 2)

    y <- matrix(0, nrow = 2, ncol = length(AT))
    for (k in 2:length(AT)) {
      y[, k] <- Ae %*% y[, k-1] + AeB * AT[k]
    }

    DT <- y[1, ]

    PSA[n] <- max(abs(DT)) * Wn^2
    PSV[n] <- max(abs(DT)) * Wn
    SD[n] <- max(abs(DT))
  }

  return(list(Tn=Tn,PSA = PSA, PSV = PSV, SD = SD))

}
