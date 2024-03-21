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
#' @param Detrend_AT boolean
#' @param Detrend_VT boolean
#' @param Detrend_DT boolean
#' @param Rebuild_AT boolean
#' @param RemoveFirstIMF_AT boolean
#' @param RemoveLastIMF_AT boolean
#' @param RemoveFirstIMF_VT boolean
#' @param RemoveLastIMF_VT boolean
#' @param RemoveFirstIMF_DT boolean
#' @param RemoveLastIMF_DT boolean
#' @param taper integer
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer

#'
#' @return data.table
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
    Detrend_AT = TRUE,
    Detrend_VT = FALSE,
    Detrend_DT = FALSE,
    PadZeros=TRUE,
    RemoveFirstIMF_AT = FALSE,
    RemoveLastIMF_AT = TRUE,
    RemoveFirstIMF_VT = FALSE,
    RemoveLastIMF_VT = TRUE,
    RemoveFirstIMF_DT = FALSE,
    RemoveLastIMF_DT = TRUE,
    Rebuild_AT = FALSE,
    TargetUnits = "mm",
    taper=2,#c(0,1,2,3) 0:None, 1:Amplitude, 2:Intensity, 3:Both
    NW = 2048,
    OVLP = 75) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)

  . <- NULL

  OK <- is.data.table(x) && !is.null(dt) && !is.null(UN)
  stopifnot(OK)

  if(taper==0) {
    .taper <- function(x){x}
    FlatZeros <- FALSE}

  if(taper==1) .taper <- function(x){.taperA(x)}
  if(taper==2) .taper <- function(x){.taperI(x)}
  if(taper==3) .taper <- function(x){.taperA(x)*.taperI(x)}

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
  if (Detrend_AT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  if (FlatZeros == TRUE) {
    Wo <- AT[,.(sapply(.SD, function(x) {.taper(x)}))]
    AT <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * Wo[[i]])]
    names(AT) <- OCID
  }






  ## Detrend AT ----
  if (Detrend_AT) {
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

  if (Detrend_AT) {
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
  # NP <-  nrow(AT)
  # ts <- seq(0,dt*(NP-1),dt)
  if(RemoveFirstIMF_AT || RemoveLastIMF_AT){
    AT <- AT[,lapply(.SD,function(x){
      # AUX <- buildIMF(t=ts,s=x,model="ceemd",trials=2)
      AUX <- buildIMF(dt=dt,s=x,model="eemd",boundary="wave",stop.rule="type5")

      nimf <- AUX$nimf
      xcols <- c(ifelse(RemoveFirstIMF_AT,1,0),ifelse(RemoveLastIMF_AT,nimf,0))
      imf <-  AUX$imf[,-xcols,with = FALSE]
      return(rowSums(imf))
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

  ## Build  Filter ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  FsNYQ <- Fs / 2 # 16
  LP <- .buildLowPassButtterworth(f = fs, Fstop = round(Fstop_LP / df) * df, Fpass = round(Fpass_LP / df) * df, Astop = 0.001, Apass = 0.95)

  HI <- .buildIntegrateFilter(f = fs) ## Integrate Filter
  HD <- .buildDerivateFilter(f = fs) ## Derivate Filter
  ## Integrate AT ----
  VT <- AT[, lapply(.SD, function(x) {
    x <- .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
  })]
  names(VT) <- OCID

  ## VT EEMD ----
  if(RemoveFirstIMF_VT || RemoveLastIMF_VT){
    VT <- VT[,lapply(.SD,function(x){
      # AUX <- buildIMF(t=ts,s=x,model="ceemd",trials=2)
      AUX <- buildIMF(dt=dt,s=x,model="eemd",boundary="wave",stop.rule="type5")

      nimf <- AUX$nimf
      xcols <- c(ifelse(RemoveFirstIMF_AT,1,0),ifelse(RemoveLastIMF_AT,nimf,0))
      imf <-  AUX$imf[,-xcols,with = FALSE]
      return(rowSums(imf))
    })]
  }
  ## Detrend VT ----
  if (Detrend_VT) {
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  ## Integrate VT ----
  DT <- VT[, lapply(.SD, function(x) {
    x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
  })]
  names(DT) <- OCID


  ## DT EEMD ----

  if(RemoveFirstIMF_DT || RemoveLastIMF_DT){
    DT <- DT[,lapply(.SD,function(x){
      # AUX <- buildIMF(t=ts,s=x,model="ceemd",trials=2)
      AUX <- buildIMF(dt=dt,s=x,model="eemd",boundary="wave",stop.rule="type5")

      nimf <- AUX$nimf
      xcols <- c(ifelse(RemoveFirstIMF_AT,1,0),ifelse(RemoveLastIMF_AT,nimf,0))
      imf <-  AUX$imf[,-xcols,with = FALSE]
      return(rowSums(imf))
    })]

  }
  ## Detrend DT ----
  if (Detrend_DT) {
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
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
  ## Rebuild AT   ----
  if(Rebuild_AT){


    VT <- DT[, lapply(.SD, function(x){
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]

    if (Detrend_VT) {
      VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    }
    names(VT) <- OCID
    AT <- VT[, lapply(.SD, function(x){
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]


  }
  if (Detrend_AT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  names(AT) <- OCID





  ## Homogenize rows ----

  NMX <- min(nrow(AT), nrow(VT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]

  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    Wo <- AT[,.(sapply(.SD, function(x) {.taper(x)}))]

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
  AUX <- data.table(ts=ts, AT = AT, VT = VT, DT = DT)
  ivars <- c("ts")
  mvars <- colnames(AUX[, -c("ts")])
  AUX <- data.table::melt(AUX, id.vars = ivars, measure.vars = mvars) |> na.omit()
  TSL <- AUX[,.(t=ts,s=value,ID=gsub("\\..*$", "", variable), OCID=gsub("^[^.]*\\.", "", variable))]
  ## Trim Records,  ----
  if(TrimZeros){
    #Warning: Different time scales from now on
    TSL <- TSL[,.trimZeros(.SD),by=c("OCID","ID")]}

  ## Return ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  return(list(TSL=TSL,Wo=Wo,PGAo=PGAo,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = UN))
}
