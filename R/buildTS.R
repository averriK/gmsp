#' buildTS
#'
#' @param x data.table
#' @param dt numeric
#' @param UN character
#' @param Fmax integer
#' @param FlatZeros boolean
#' @param PadZeros boolean
#' @param TrimZeros boolean
#' @param Resample boolean
#' @param Detrend.AT boolean
#' @param Detrend.VT boolean
#' @param Detrend.DT boolean
#' @param LowPass.AT boolean
#' @param LowPass.VT boolean
#' @param LowPass.DT boolean
#' @param Rebuild boolean
#' @param EMD.AT boolean
#' @param EMD.VT boolean
#' @param EMD.DT boolean
#' @param EMD.method string
#' @param removeIMF1.AT boolean
#' @param removeIMFn.AT boolean
#' @param removeIMF1.VT boolean
#' @param removeIMFn.VT boolean
#' @param removeIMF1.DT boolean
#' @param removeIMFn.DT boolean
#' @param Taper integer
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
    x, dt, UN,
    Fmax = 25,
    Resample = TRUE,
    FlatZeros = TRUE,
    TrimZeros = TRUE,
    Detrend.AT = FALSE,
    Detrend.VT = FALSE,
    Detrend.DT = FALSE,
    LowPass.AT = FALSE,
    LowPass.VT = FALSE,
    LowPass.DT = FALSE,
    PadZeros=TRUE,
    EMD.AT = TRUE,
    EMD.VT = TRUE,
    EMD.DT = TRUE,
    EMD.method ="emd",
    removeIMF1.AT = 0,
    removeIMFn.AT = 0,
    removeIMF1.VT = 0,
    removeIMFn.VT = 0,
    removeIMF1.DT = 0,
    removeIMFn.DT = 0,
    Rebuild = FALSE,
    TargetUnits = "mm",
    Taper=1,#c(0,1,2,3) 0:None, 1:Amplitude, 2:Intensity, 3:Both
    NW = 2048,
    OVLP = 75) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)

  . <- NULL

  OK <- is.data.table(x) && !is.null(dt) && !is.null(UN)
  stopifnot(OK)

  if(Taper==0) {
    .taper <- function(x){x}
    FlatZeros <- FALSE}

  if(Taper==1) .taper <- function(x){.taperA(x)}
  if(Taper==2) .taper <- function(x){.taperI(x)}
  if(Taper==3) .taper <- function(x){.taperA(x)*.taperI(x)}

  ATo <- copy(x)
  ## Check Record ----
  NP <- nrow(ATo)
  if (  NP == 0 ||dt == 0 ||any(is.na(ATo)) ||max(abs(ATo)) == 0 ) return(NULL)
  if (NP < NW) return(NULL)
  OCID <- names(ATo)
  if (length(unique(OCID)) < 3) return(NULL)


  ## Scale Units  ----

  if (grepl(UN, pattern = "[///+]")) {
    UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
  }
  if (!(tolower(UN) %in% c("mm", "cm", "m", "gal", "g"))) return(NULL)

  if (tolower(UN) != TargetUnits) {
    SFU <- .getSF(SourceUnits = tolower(UN), TargetUnits = TargetUnits)
    # AT <- map(AT,function(x){x*SFU})
    # ATo[, (colnames(ATo)) := lapply(.SD, function(x) {x * SFU})]
    ATo <- ATo[,.(sapply(.SD, function(x) {x * SFU}))]
  }
  Fs <- 1 / dt




  ## Set Scale Reference ----
  DUMMY <- NULL
  PGAo <- apply(ATo, 2, function(x) { max(abs(x))})
  ## Scale record ----
  AT <-ATo[, .(sapply(.SD, function(x){x/max(abs(x))}))]



  ## Flat Zeros (A+I)  ----
  if (Detrend.AT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  if (FlatZeros == TRUE) {
    Wo <- AT[,.(sapply(.SD, function(x) {.taper(x)}))]
    AT <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * Wo[[i]])]
    names(AT) <- OCID
  }






  ## Detrend AT ----
  if (Detrend.AT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
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

    AT <- AT[, lapply(.SD, function(x) {
      x <- ffilter(wave = x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    AT <- AT[, lapply(.SD, function(x) {
      x <- signal::resample(x, TargetFs, Fs)
      return(x)
    })]
    Fs <- TargetFs
    names(AT) <- OCID
  }
  ## Detrend AT ----
  if (Detrend.AT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    # names(AT) <- OCID
  }
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  dt <- 1/Fs
  NP <- nrow(AT)

  ## Flat Zeros (A+I)  ----

  if (FlatZeros == TRUE) {
    Wo <- AT[,.(sapply(.SD, function(x) {.taper(x)}))]
    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    names(AT) <- OCID
  }



  ## AT EEMD ----
  if (EMD.AT) {
    AT <- AT[,lapply(.SD,function(x){
      AUX <- buildIMF(dt=dt,s=x,method=EMD.method,max.imf=10)
      nimf <- AUX$nimf
      i <- removeIMF1.AT
      j <- removeIMFn.AT
      COLS <- colnames(AUX$imf)[(i+1):(nimf-j)]
      x <-  AUX$imf[,COLS,with = FALSE] |> rowSums()
      return(x)
    })]
  }


  ## Padding Zeros ----
  if(PadZeros){
    NP <- nrow(AT)
    NZ <- .getNZ(NP)
    if (NZ > 0) {
      O <- data.table()[, (colnames(AT)) := list(rep(0, NZ))]
    } else {
      O <- data.table()[, (colnames(AT)) := list(rep(0, NW))]
    }
    AT <- rbindlist(list(O, AT))

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
  ## AT-> VT ----
  VT <- AT[, lapply(.SD, function(x) {
    x <- .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW
  })]
  names(VT) <- OCID
  ## Lowpass VT ----
  if(LowPass.VT){
    VT <- VT[, lapply(.SD, function(x) {
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    names(VT) <- OCID
  }


  ## EMD VT ----
  if (EMD.VT) {
    VT <- VT[,lapply(.SD,function(x){
      # AUX <- buildIMF(t=ts,s=x,eemd.lib="ceemd",trials=2)
      AUX <- buildIMF(dt=dt,s=x,method=EMD.method,max.imf=15)

      nimf <- AUX$nimf
      i <- removeIMF1.VT
      j <- removeIMFn.VT
      COLS <- colnames(AUX$imf)[(i+1):(nimf-j)]
      x <-  AUX$imf[,COLS,with = FALSE] |> rowSums()
      return(x)
    })]
  }
  ## Detrend VT ----
  if (Detrend.VT) {
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  ## Integrate VT ----
  DT <- VT[, lapply(.SD, function(x) {
    x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)
  })]
  names(DT) <- OCID


  ## Lowpass DT ----
  if(LowPass.DT){
    DT <- DT[, lapply(.SD, function(x) {
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    names(DT) <- OCID
  }

  ## EMD DT ----
  if (EMD.DT) {
    DT <- DT[,lapply(.SD,function(x){
      AUX <- buildIMF(dt=dt,s=x,method=EMD.method,max.imf=15)

      nimf <- AUX$nimf
      i <- removeIMF1.DT
      j <- removeIMFn.DT
      COLS <- colnames(AUX$imf)[(i+1):(nimf-j)]
      x <-  AUX$imf[,COLS,with = FALSE] |> rowSums()
      return(x)
    })]
  }
  ## Detrend DT ----
  if (Detrend.DT) {
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }



  ## Rebuild ----
  if(Rebuild){
    # browser()
    # Flat Zeros on DT based on AT
    NMX <- min(nrow(AT), nrow(DT))
    AT <- AT[-((NMX):.N)]
    DT <- DT[-((NMX):.N)]
    Wo <- AT[,.(sapply(.SD, function(x) {.taperA(x)*.taperI(x)}))]
    idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
    DT <- DT[idx]
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
    ## Derivate DT ----
    VT <- DT[, lapply(.SD, function(x){
      .derivate(s=x,dt=dt) })]
    if (Detrend.VT) {
      VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    }

    ## Derivate VT ----
    AT <-  VT[, lapply(.SD, function(x){
      .derivate(s=x,dt=dt) })]

    if (Detrend.AT) {
      AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    }
    NMX <- min(nrow(AT), nrow(DT))
    AT <- AT[-((NMX):.N)]
    VT <- VT[-((NMX):.N)]
    DT <- DT[-((NMX):.N)]

  }
  ## Homogenize rows ----

  NMX <- min(nrow(AT), nrow(VT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]
  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    # browser()
    Wo <- AT[,.(sapply(.SD, function(x) {.taperA(x)}))]

    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]

    names(AT) <- OCID
    names(VT) <- OCID
    names(DT) <- OCID
  }




  ## Restore Scale ----

  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
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

  return(list(TSL=TSL,TSW=TSW,Wo=Wo,PGAo=PGAo,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = UN))
}
