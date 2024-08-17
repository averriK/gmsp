#' Title
#'
#' @param .x data.table
#' @param s vector
#' @param t vector
#' @param dt double
#' @param Tn vector
#' @param xi double
#'
#' @return list
#' @export 
#'
#' @import expm
#' @import data.table
#' @importFrom pracma pinv
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
  # if(is.null(Tn)){
  #   seq1 <- (10^seq(log10(0.001), log10(0.1), length.out = 40)) |> round(3)
  #   seq2 <- (10^seq(log10(0.1), log10(1), length.out = 20)) |> round(2)
  #   seq3 <- (10^seq(log10(1), log10(10), length.out = 10)) |> round(1)
  #   Tn <- c(0,seq1, seq2, seq3) |> unique() 
  # }
  
  if(is.null(Tn)){
    Tn <- c(0.01, 0.02, 0.022, 0.025, 0.029, 0.03, 0.032, 0.035, 0.036,
            0.04, 0.042, 0.044, 0.045, 0.046, 0.048, 0.05, 0.055, 0.06, 0.065, 0.067, 0.07,0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.11, 0.12, 0.125, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 0.3, 0.32, 0.34,0.35, 0.36, 0.38, 0.4, 0.42, 0.44, 0.45, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 0.667,0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,2, 2.2, 2.4, 2.5, 2.6, 2.8, 3, 3.2, 3.4, 3.5, 3.6, 3.8, 4, 4.2, 4.4, 4.6, 4.8, 5, 7.5, 10)
  }
  
  
  PGA <- max(abs(AT))
  PSA <- numeric(length(Tn))
  for (i in seq_along(Tn)) {
    omega_n <- 2 * pi / Tn[i]
    C <- 2 * xi * omega_n
    K <- omega_n^2
    y <- matrix(0, nrow = 2, ncol = length(AT))
    A <- matrix(c(0, 1, -K, -C), nrow = 2, byrow = TRUE)  # Matrix definition
    Ae <- expm(A * dt)  # Matrix exponentiation
    AeB <- (Ae - diag(2)) %*% pracma::pinv(A) %*% matrix(c(0, 1), nrow = 2,byrow=FALSE)  # Compute AeB
    
    for (k in 2:length(AT)) {
      y[, k] <- Ae %*% y[, k-1] + AeB * AT[k]
    }
    
    displ_max <- max(abs(y[1, ]))  # Spectral relative displacement (m)
    PSA[i] <- displ_max * K #omega_n^2  # Pseudo spectral acceleration (m/s2)
  }
  Tn <- c(0,Tn)
  PSA <- c(PGA,PSA)
  ## Pack Spectra  ----
  
  return(list( Tn=Tn, xi=xi,PSA = PSA))
  
}

# seq1 <- (10^seq(log10(0.001), log10(0.1), length.out = 40)) |> round(3)
# seq2 <- (10^seq(log10(0.1), log10(1), length.out = 20)) |> round(2)
# seq3 <- (10^seq(log10(1), log10(10), length.out = 10)) |> round(1)
# c(seq1, seq2, seq3) |> unique() |> print()
