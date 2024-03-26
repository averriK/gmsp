#' buildTS
#'
#' @param x data.table
#' @param dt numeric
#' @param ts numeric
#' @param Order integer
#' @param UN character
#' @param Fmax integer
#' @param Resample boolean
#' @param LowPass boolean
#' @param removeIMF1 boolean
#' @param removeIMFn boolean
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
    LowPass = FALSE,
    removeIMF1 = 0,
    removeIMFn = 0,
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





  ## Resample ----
  if(Resample){
    TargetFs <- as.integer(5*Fmax) # 80/200 Hz
    X <- .resample(X,dt=dt,TargetFs=TargetFs,Fmax=Fmax,NW=NW,OVLP=OVLP)
    Fs <- TargetFs
    dt <- 1/Fs
  }


  # Case #2. Acceleration Time Histories
  if(Order==2){
    AT <- X
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(AT) <- OCID

    VT <- AT[, lapply(.SD, function(x){ .integrate(dX=x,dt=dt) })]
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(VT) <- OCID

    DT <- VT[, lapply(.SD, function(x){ .integrate(dX=x,dt=dt) })]
    # DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(DT) <- OCID
  }

  # Case #1. Velocity Time Histories
  if(Order==1){
    VT <- X
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]

    DT <- VT[, lapply(.SD, function(x){.integrate(dX=x,dt=dt) })]
    # DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(DT) <- OCID


  }

  # Case #0. Displacement Time Historyes
  if(Order==0){
    DT <- X


  }
  NMX <- min(nrow(AT),nrow(VT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]
  ## REBUILD ----
  Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
  # Wo <- DT[,.(sapply(.SD, function(x) {.taperI(x)}))]

  # idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
  # DT <- DT[idx]
  DT <-DT[, .(sapply(.SD, function(x){.detrend(X=x,dt=dt,removeIMF1=removeIMF1,removeIMFn=removeIMFn)}))]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]

  ## Derivate DT ----
  VT <- DT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
  VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]

  ## Derivate VT ----
  AT <-  VT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
  AT <- AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  NMX <- min(nrow(AT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]
  # # .taperI muy restrictivo para aceleraciones
  # Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
  Wo <- AT[,.(sapply(.SD, function(x) {.taperA(x,Astop=1e-5,Apass=1e-4)}))]
  idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
  AT <- AT[idx]
  VT <- VT[idx]
  DT <- DT[idx]
  # AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  # VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  # DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]


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
