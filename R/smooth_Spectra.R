#' Title
#'
#' @param Tn vector. Frequency (Periods) vector
#' @param A vector. Amplitude spectrum
#' @param method string
#' @param window numeric
#' @param po numeric
#' @param cv boolean
#' @param sigma numeric
#'
#' @return list
#' @export 
#'
#' @import data.table
#' @importFrom stats filter
#' @importFrom signal sgolayfilt
#' @importFrom zoo rollmean
#' @importFrom stats smooth.spline
#' @importFrom stats predict
#' @importFrom stats dnorm
#' @importFrom stats spline
#'
#' 
smooth_Spectra <- function(Tn,A,method="none",window=5,po=3,cv=FALSE,sigma=2){
  stopifnot(tolower(method) %in% c("none","ma","sg","sg","gk","ema","sm"))
  As <- switch(
    method,
    "ma" = {
      stats::filter(A, rep(1/window, window), sides=2)
    },
    "sg" = {
      signal::sgolayfilt(A, p=po, n=ifelse(window %% 2 == 1, window, window + 1))
    },
    "gk" = {
      kernel <- stats::dnorm(seq(-window, window, length.out=2*window+1), sd=sigma)
      kernel <- kernel / sum(kernel)  # Normalize kernel
      stats::convolve(A, kernel, type='filter')
    },
    "ema" = {
      zoo::rollmean(A, window, align='center', fill=NA)  # Exponential Moving Average
    },
    "sm" = {
      spline_fit <- stats::smooth.spline(Tn, A,cv=cv)
      predict(spline_fit, Tn)$y
    },
    "none" = {
      A
    }
  )

  # Ensure the smoothed amplitude spectrum has the same length as the frequency vector
  As <- c(As, rep(0, length(Tn) - length(As)))
  DT <- list(As=As, Tn=Tn)
}
