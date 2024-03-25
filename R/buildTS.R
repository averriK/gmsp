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
#' @param EMD.AT boolean
#' @param EMD.VT boolean
#' @param EMD.DT boolean
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
    Fmax = 25,
    Resample = TRUE,
    FlatZeros = TRUE,
    TrimZeros = TRUE,
    Detrend = FALSE,
    LowPass = FALSE,
    EMD.AT = FALSE,
    EMD.VT = FALSE,
    EMD.DT = FALSE,
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
  X <- copy(x)

  NP <- ifelse(is.data.table(X),nrow(X),ifelse(is.vector(X),length(X)))
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


  # ATo <- copy(x)
  OCID <- names(X)




  ## Set Scale Reference ----
  DUMMY <- NULL
  Ao <- apply(X, 2, function(x) { max(abs(x))})
  ## Scale record ----
  X <-X[, .(sapply(.SD, function(x){x/max(abs(x))}))]
  X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]


  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    Wo <- X[,.(sapply(.SD, function(x) {.taperI(x)}))]
    X <- X[, lapply(seq_along(.SD), function(i) .SD[[i]] * Wo[[i]])]
    names(X) <- OCID
    X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]
  }


  ## Resample ----
  if(Resample){
    TargetFs <- as.integer(4*Fmax) # 160/200 Hz
    Fpass_LP <- Fmax # 20/25 Hz
    Fstop_LP <- 1.2*Fmax # 25/30 Hz
    Fs <- 1/dt #
    df <- Fs / NW #
    fs <- seq(from = 0, by = df, length.out = NW / 2)
    LP <- .buildLowPassButtterworth(f = fs, Fstop = round(1 * Fstop_LP / df) * df, Fpass = round(1 * Fpass_LP / df) * df, Astop = 0.001, Apass = 0.95)


    X <- X[, lapply(.SD, function(x) {
      ffilter(wave = x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
    })]
    X <- X[, lapply(.SD, function(x) {
      signal::resample(x, TargetFs, Fs)

    })]
    Fs <- TargetFs
    names(X) <- OCID
    X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]
  }




  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  dt <- 1/Fs
  NP <- nrow(X)

  ## Flat Zeros (A+I)  ----

  if (FlatZeros == TRUE) {
    Wo <- X[,.(sapply(.SD, function(x) {.taperI(x)}))]
    X <- X[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    names(X) <- OCID
  }







  ## Build  Filters ----
  Fs <- 1/dt #
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  Fpass_LP <- Fmax # 20/25 Hz
  Fstop_LP <- 1.2*Fmax # 25/30 Hz
  LP <- .buildLowPassButtterworth(f = fs, Fstop = round(Fstop_LP / df) * df, Fpass = round(Fpass_LP / df) * df, Astop = 0.01, Apass = 0.99)


  HI <- .buildIntegrateFilter(f = fs) ## Integrate Filter
  HD <- .buildDerivateFilter(f = fs) ## Derivate Filter
  X <- .padZeros(X)
  X <- X[, lapply(.SD, function(x) {
    seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
  })]
  X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]
  names(X) <- OCID



  ## AT-> VT ----
  # if X is AT, integrate velocities and displacements
  if(Order==2){ # Accelerations
    AT <- X
    VT <- AT[, lapply(.SD, function(x) {
      .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW
    })]
    VT <- VT[, lapply(.SD, function(x) {
      seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)})]
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(VT) <- OCID
    DT <- VT[, lapply(.SD, function(x) {
      .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW
    })]

  }

  if(Order==1){# Velocities
    VT <- X

    Wo <- VT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
    VT <- VT[idx]

    DT <- VT[, lapply(.SD, function(x) {
      .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW
    })]
    DT <- DT[, lapply(.SD, function(x) {
      seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)})]
    names(DT) <- OCID
  }

  if(Order==0){# Displacements
    DT <- X
  }
  # Wo <- AT[,.(sapply(.SD, function(x) {.taperA(x)*.taperI(x)}))]
  # idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
  # AT <- AT[idx]

  ## EMD DT ----
  if (EMD.DT) {
    DT <- DT[,lapply(.SD,function(x){
      AUX <- buildIMF(dt=dt,s=x,method="emd",max.imf=15)

      nimf <- AUX$nimf
      i <- removeIMF1.DT
      j <- removeIMFn.DT
      COLS <- colnames(AUX$imf)[(i+1):(nimf-j)]
      x <-  AUX$imf[,COLS,with = FALSE] |> rowSums()
      return(x)
    })]
  }


  ## Rebuild ----
  if(!Rebuild){
    NMX <- min(nrow(AT), nrow(DT))
    AT <- AT[-((NMX):.N)]
    VT <- VT[-((NMX):.N)]
    DT <- DT[-((NMX):.N)]

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

    # NMX <- min(nrow(AT), nrow(DT))
    # AT <- AT[-((NMX):.N)]
    # DT <- DT[-((NMX):.N)]
    Wo <- DT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
    DT <- DT[idx]
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
    ## Derivate DT ----
    VT <- DT[, lapply(.SD, function(x){
      .derivate(s=x,dt=dt) })]
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]

    ## Derivate VT ----
    AT <-  VT[, lapply(.SD, function(x){
      .derivate(s=x,dt=dt) })]
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  ## Homogenize rows ----

  # NMX <- min(nrow(AT), nrow(VT), nrow(DT))
  # AT <- AT[-((NMX):.N)]
  # VT <- VT[-((NMX):.N)]
  # DT <- DT[-((NMX):.N)]
  ## Flat Zeros (A+I)  ----

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
  if(TrimZeros){
    #Warning: Different time scales from now on
    TSL <- TSL[,.trimZeros(.SD),by=c("OCID","ID")]}

  ## Return ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  return(list(TSL=TSL,TSW=TSW,Wo=Wo,Ao=Ao,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = UN))
}
