
#' Title
#' 
#' @param .x data.table
#' @param s numeric vector
#' @param ts numeric vector
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
get_imf <- function(.x=NULL, s=NULL, ts=NULL, dt=NULL, method="emd", boundary="wave", max.imf=15, noise.type="gaussian", noise.amp=0.5e-7, trials=10, stop.rule="type5", verbose=FALSE) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)
  
  stopifnot(tolower(method) %in% c("emd", "eemd", "ceemd"))
  stopifnot(tolower(stop.rule) %in% c("type1", "type2", "type3", "type4", "type5"))
  stopifnot(tolower(boundary) %in% c("none", "wave", "symmetric", "periodic", "evenodd"))
  stopifnot(tolower(noise.type) %in% c("uniform", "gaussian"))

  if (is.null(s) && is.null(.x)) {
    stop("No acceleration data provided.")
  }
  if (!is.null(.x)) {
    s <- .x$s
    ts <- .x$t
  }
  if (is.null(ts) && is.null(.x) && is.null(dt)) {
    stop("No time data provided.")
  }
  
  if (!is.null(ts)) {
    dt <- mean(diff(ts))
    to <- min(ts)
  }
  
  ts <- seq(0, (length(s)-1)*dt, by=dt)
  
  DT <- data.table(t=ts, s=s)
  DT <- .trimZeros(DT)
  n <- nrow(DT)
  rm(s)
  
  if (tolower(method) == "emd") {
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf, stoprule=stop.rule)
    IMF <- as.data.table(AUX$imf)
    RES <- AUX$residue
  }
  
  if (tolower(method) == "eemd") {
    nimf <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf, stoprule=stop.rule)$nimf
    DIR <- tempdir(check = TRUE)
    hht::EEMD(sig=DT$s, tt=DT$t, nimf=nimf, max.imf=max.imf, boundary=boundary, noise.amp=noise.amp, noise.type=noise.type, trials=trials, stop.rule=stop.rule, trials.dir=DIR, verbose=verbose)
    AUX <- EEMDCompile(trials.dir=DIR, trials=trials, nimf=nimf) |> suppressWarnings()
    unlink(DIR, force=TRUE, recursive=TRUE)
    IMF <- as.data.table(AUX$averaged.imfs)
    RES <- AUX$averaged.residue
  }
  
  if (tolower(method) == "ceemd") {
    AUX <- hht::CEEMD(sig=DT$s, tt=DT$t, noise.amp=noise.amp, noise.type=noise.type, trials=trials, stop.rule=stop.rule, verbose=verbose)
    IMF <- as.data.table(AUX$imf)
    RES <- AUX$residue
  }
  
  if(ncol(IMF) == 0) {
    warning("No IMFs were found.")
    return(NULL)
  }
  # Zero padding to ensure IMF has the same length as s and t
  if (nrow(IMF) < nrow(DT)) {
    pad_length <- nrow(DT) - nrow(IMF)
    padding <- data.table(matrix(0, nrow=pad_length, ncol=ncol(IMF)))
    setnames(padding, names(IMF))
    IMF <- rbind(IMF, padding)
  }
  
  # Rename columns before adding signal and residue
  setnames(IMF, c(paste0("IMF", seq_len(ncol(IMF)))))
  
  # Add signal and residue as the last columns of IMF
  TSW <- data.table(t=DT$t+to, signal=DT$s, IMF,residue=RES)
  #
  ivars <- c("t")
  mvars <- colnames(TSW[, -c("t")])
  AUX <- melt(TSW, id.vars=ivars, measure.vars=mvars)
  TSL <- AUX[, list(t=t, s=value, IMF=variable)]
  
  return(TSL)
}
