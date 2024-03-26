#' buildTS
#'
#' @param x data.table
#' @param dt numeric
#' @param ts numeric
#' @param Order integer
#' @param UN character
#' @param Fmax integer
#' @param FlatZeros boolean
#' @param TrimZeros boolean
#' @param Resample boolean
#' @param Detrend boolean
#' @param LowPass boolean
#' @param Rebuild boolean
#' @param removeIMF1.AT boolean
#' @param removeIMFn.AT boolean
#' @param removeIMF1.VT boolean
#' @param removeIMFn.VT boolean
#' @param removeIMF1.DT boolean
#' @param removeIMFn.DT boolean
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer

#'
#' @return list
#' @export buildTS
#'
#' @examples
#'
#' @importFrom data.table :=
#' @importFrom stats na.omit
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom seewave ffilter
#' @importFrom signal resample
#' @importFrom stringr str_split
#' @importFrom purrr map
#'
#'
buildTS <- function(
    x, ts=NULL ,dt=NULL, UN,
    Order=2,
    Fmax = 16,
    Resample = TRUE,
    FlatZeros = TRUE,
    TrimZeros = TRUE,
    Detrend = FALSE,
    LowPass = FALSE,
    removeIMF1.AT = 0,
    removeIMFn.AT = 0,
    removeIMF1.VT = 0,
    removeIMFn.VT = 0,
    removeIMF1.DT = 0,
    removeIMFn.DT = 0,
    Rebuild = FALSE,
    TargetUnits = "mm",
    NW = 1024,
    OVLP = 75) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)

  . <- NULL
  X <- copy(x) |> as.data.table()

  NP <- nrow(X)
  stopifnot(NP>4 && NP>=NW)
  if(!is.null(ts)){
    dt <- diff(ts) |> mean()
  }
  ts <- seq(0,(NP-1)*dt,by=dt)
  Fs <- 1 / dt

  ## Scale Units  ----

  if (grepl(UN, pattern = "[///+]")) {
    UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
  }
  if (!(tolower(UN) %in% c("mm", "cm", "m", "gal", "g"))) return(NULL)

  if (tolower(UN) != TargetUnits) {
    SFU <- .getSF(SourceUnits = tolower(UN), TargetUnits = TargetUnits)
    # AT <- map(AT,function(x){x*SFU})
    # ATo[, (colnames(ATo)) := lapply(.SD, function(x) {x * SFU})]
    X <- X[,.(sapply(.SD, function(x) {x * SFU}))]
  }

  OCID <- names(X)




  ## Set Scale Reference ----
  Ao <- apply(X, 2, function(x) { max(abs(x))})
  ## Scale record ----
  X <-X[, .(sapply(.SD, function(x){x/max(abs(x))}))]
  # Detrend
  X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]


  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    Wo <- X[,.(sapply(.SD, function(x) {.taperI(x)}))]
    X <- X[, lapply(seq_along(.SD), function(i) .SD[[i]] * Wo[[i]])]
    names(X) <- OCID

  }


  ## Resample ----
  if(Resample){
    TargetFs <- as.integer(5*Fmax) # 80/200 Hz
    X <- .resample(X,Fs=Fs,Fmax=Fmax,NW=NW,OVLP=OVLP)
    Fs <- TargetFs
    dt <- 1/Fs
  }

#
#
#
#   df <- Fs / NW # 0.03125#
#   fs <- seq(from = 0, by = df, length.out = NW / 2)
#   dt <- 1/Fs
#   NP <- nrow(X)

  ## Flat Zeros (A+I)  ----

  if (FlatZeros == TRUE) {
    Wo <- X[,.(sapply(.SD, function(x) {.taperI(x)}))]
    X <- X[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    names(X) <- OCID
  }

  # Trim Zeros
  Wo <- X[,.(sapply(.SD, function(x) {.taperI(x)}))]
  idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
  X <- X[idx]

  # Case #1. Acceleration Time Histories
  if(Order==2){
    AT <- X
    VT <- AT[, lapply(.SD, function(x){.integrate(dX=x,dt=dt) })]
    DT <- VT[, lapply(.SD, function(x){.integrate(dX=x,dt=dt) })]
  }

  # Case #2. Velocity Time Histories
  if(Order==1){
    VT <- X
    DT <- VT[, lapply(.SD, function(x){.integrate(dX=x,dt=dt) })]
  }

  # Case #3. Displacement Time Historyes
  if(Order==0){
    DT <- X
  }
  NMX <- min(nrow(AT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]



  ## Rebuild ----
  if(!Rebuild){

    Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
    AT <- AT[idx]
    VT <- VT[idx]
    DT <- DT[idx]

    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  if(Rebuild){
    Wo <- DT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
    DT <- DT[idx]
    # EMD DT ----

    ## Derivate DT ----
    VT <- DT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
    # EMD cT ----

    ## Derivate VT ----
    AT <-  VT[, lapply(.SD, function(x){.derivate(s=x,dt=dt) })]


  }

  ## Restore Scale ----

  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
  names(AT) <- OCID
  names(VT) <- OCID
  names(DT) <- OCID

  ## Pack Time Series  ----
  NP <-  nrow(AT)
  ts <- seq(0,dt*(NP-1),dt)
  TSW <- data.table(ts=ts, AT = AT, VT = VT, DT = DT)
  ivars <- c("ts")
  mvars <- colnames(TSW[, -c("ts")])
  AUX <- data.table::melt(TSW, id.vars = ivars, measure.vars = mvars) |> na.omit()
  TSL <- AUX[,.(t=ts,s=value,ID=gsub("\\..*$", "", variable), OCID=gsub("^[^.]*\\.", "", variable))]
  ## Trim Records,  ----
  # if(TrimZeros){
  #   #Warning: Different time scales from now on
  #   TSL <- TSL[,.trimZeros(.SD),by=c("OCID","ID")]}

  ## Return ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  return(list(TSL=TSL,TSW=TSW,Wo=Wo,Ao=Ao,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = UN))
}
