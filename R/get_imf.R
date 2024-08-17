
#' Title
#' 
#' @param .x data.table
#' @param s numeric vector
#' @param t numeric vector
#' @param dt numeric
#' @param method string
#' @param boundary string
#' @param max.imf integer
#' @param noise.type string
#' @param noise.amp numeric
#' @param trials integer
#' @param stop.rule string
#' @param verbose boolean
#'
#' @importFrom data.table :=
#' @importFrom stats na.omit
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom hht EEMDCompile
#' @importFrom hht CEEMD
#' @importFrom hht EEMD
#' @importFrom EMD emd
#'
#' @return data.table
#' @export 
#'
#' @examples
#'
get_imf <- function(.x=NULL, s=NULL, t=NULL, dt=NULL, method="emd", boundary="wave", max.imf=15, noise.type="gaussian", noise.amp=0.5e-7, trials=10, stop.rule="type5", verbose=FALSE) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)
  . <- NULL
  stopifnot(tolower(method) %in% c("emd", "eemd", "ceemd"))
  stopifnot(tolower(stop.rule) %in% c("type1", "type2", "type3", "type4", "type5"))
  stopifnot(tolower(boundary) %in% c("none", "wave", "symmetric", "periodic", "evenodd"))
  stopifnot(tolower(noise.type) %in% c("uniform", "gaussian"))
  
  if (is.null(s) && is.null(.x)) {
    stop("No acceleration data provided.")
  }
  if (!is.null(.x)) {
    s <- .x$s
    t <- .x$t
  }
  if (is.null(t) && is.null(.x) && is.null(dt)) {
    stop("No time data provided.")
  }
  
  if (!is.null(t)) {
    dt <- diff(t) |> mean()
  }
  
  n <- length(s)
  stopifnot(n > 4)
  
  t <- seq(0, (n-1)*dt, by=dt)
  
  DT <- data.table(t=t, s=s)
  DT <- .trimZeros(DT)
  
  if (tolower(method) == "emd") {
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf, stoprule=stop.rule)
    IMF <- AUX$imf |> as.data.table()
    RES <- AUX$residue |> as.vector() |> unname()
  }
  
  if (tolower(method) == "eemd") {
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf, stoprule=stop.rule)
    nimf <- AUX$nimf
    DIR <- tempdir(check = TRUE)
    hht::EEMD(sig=DT$s, tt=DT$t, nimf=nimf, max.imf=max.imf, boundary=boundary, noise.amp=noise.amp, noise.type=noise.type, trials=trials, stop.rule=stop.rule, trials.dir=DIR, verbose=verbose)
    AUX <- EEMDCompile(trials.dir=DIR, trials=trials, nimf=nimf) |> suppressWarnings()
    unlink(DIR, force=TRUE, recursive=TRUE)
    IMF <- AUX$averaged.imfs |> as.data.table()
    RES <- AUX$averaged.residue |> unname()
  }
  
  if (tolower(method) == "ceemd") {
    AUX <- hht::CEEMD(sig=DT$s, tt=DT$t, noise.amp=noise.amp, noise.type=noise.type, trials=trials, stop.rule=stop.rule, verbose=verbose)
    IMF <- AUX$imf |> as.data.table()
    RES <- AUX$residue |> unname()
  }
  
  # Add signal and residue as the last columns of IMF
  IMF <- cbind(IMF, signal=DT$s, residue=RES)
  # Zero padding
  i <- as.integer(DT$t[1] / dt)
  j <- as.integer(last(DT$t) / dt)
  
  OB <- data.table(matrix(0, nrow=n-j, ncol=ncol(IMF)))
  OA <- data.table(matrix(0, nrow=i-1, ncol=ncol(IMF)))
  IMF <- rbindlist(list(OA, IMF, OB), use.names=FALSE)
  
  # Set names
  names(IMF) <- c(paste0("IMF", seq_len(ncol(IMF)-2)), "signal", "residue")
  
  TSW <- data.table(t, IMF)
  ivars <- c("t")
  mvars <- colnames(TSW[, -c("t")])
  AUX <- data.table::melt(TSW, id.vars=ivars, measure.vars=mvars) |> na.omit()
  TSL <- AUX[, .(t=t, s=value, IMF=variable)]
  
  return(TSL)
}
