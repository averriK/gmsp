
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
#' @param plot boolean
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
#' @return list
#' @export 
#'
#' @examples
#'
build_EMD_plot <- function(....) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)
  
  # Call get_imf to get the IMF data
  TSL <- get_imf(...)
  
  IMF <- dcast(TSL, t ~ IMF, value.var="s")
  RES <- IMF$residue
  SIGNAL <- IMF$s
  IMF[, c("s", "residue") := NULL]  # Remove signal and residue columns to avoid duplication
  
  Tm <- IMF[, lapply(.SD, function(x) {.getTm(x, Fs=1/dt)})]
  wm <- 2*pi/Tm
  fm <- 1/Tm
  PGA <- IMF[, lapply(.SD, function(x) {max(abs(x))})]
  nimf <- ncol(IMF)
  
  # Offset for plotting
  M <- IMF[, lapply(.SD, function(col) col + seq(0, (ncol(IMF)-1) * 1.25 * ceiling(max(IMF) - min(IMF)), length.out = nrow(IMF)))]
  M <- cbind(signal=SIGNAL, M, residue=RES)
  
  AUX <- cbind(IMF[, .(t)], M)
  ivars <- c("t")
  mvars <- colnames(AUX[, -c("t")])
  DATA <- melt(AUX, id.vars=ivars, measure.vars=mvars) |> na.omit()
  DATA <- DATA[, .(X=t, Y=value, ID=variable)]
  rm(M)
  
  return(list(plot.data=DATA, plot.offset=offset))
}


