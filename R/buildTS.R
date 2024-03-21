#' buildTS
#'
#' @param x data.table
#' @param dt numeric
#' @param UN character
#' @param Fmax integer
#' @param FlatZeros boolean
#' @param PadZeros boolean
#' @param TrimZeros boolean
#' @param Modes numeric
#' @param Resample boolean
#' @param DetrendAT boolean
#' @param DetrendVT boolean
#' @param DetrendDT boolean
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
    x, dt, UN, Modes=NULL,
    Fmax = 25,
    Resample = TRUE,
    FlatZeros = TRUE,
    TrimZeros = TRUE,
    DetrendAT = TRUE,
    DetrendVT = TRUE,
    DetrendDT = TRUE,
    PadZeros=TRUE,
    TargetUnits = "mm",
    NW = 2048,
    OVLP = 75) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)

  . <- NULL

  OK <- is.data.table(x) && !is.null(dt) && !is.null(UN)
  stopifnot(OK)

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
  # browser()
  if (FlatZeros == TRUE) {
    # WoA <- AT[,.(sapply(.SD, function(x) {.taperA(x)}))]
    # WoI <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)*.taperA(x)}))]
    AT <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * Wo[[i]])]
    names(AT) <- OCID
  }
  # names(AT) <- OCID
  if (DetrendAT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    # names(AT) <- OCID
  }



  ## Resample ----
  if(Resample){
    TargetFsNYQ <- as.integer(4*Fmax) # 80/100 Hz
    TargetFs <- 2*TargetFsNYQ # 160/200 Hz
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
  # browser()
  if (DetrendAT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    # names(AT) <- OCID
  }
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  dt <- 1/Fs
  NP <- nrow(AT)

  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    # WoA <- AT[,.(sapply(.SD, function(x) {.taperA(x)}))]
    # WoI <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)*.taperA(x)}))]
    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    names(AT) <- OCID
  }




  ## Padding Zeros ----
  if(PadZeros){
    NP <- nrow(AT)
    NZ <- .getNZ(NP)
    if (NZ > 0) {
      ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NZ))]
    } else {
      ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NW))]
    }
    AT <- rbindlist(list(ZEROS, AT))

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
  if (DetrendVT) {
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }

  ## Integrate VT ----
  DT <- VT[, lapply(.SD, function(x) {
    x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
  })]
  names(DT) <- OCID
  if (DetrendDT) {
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }




  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    # WoA <- AT[,.(sapply(.SD, function(x) {.taperA(x)}))]
    Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)*.taperA(x)}))]
    AT <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * Wo[[i]])]
    names(AT) <- OCID
  }


  ## Homogeinize rows ----
  # browser()
  NMX <- min(nrow(AT), nrow(VT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]

  ## Flat Zeros (A+I)  ----
  if (FlatZeros == TRUE) {
    # WoA <- AT[,.(sapply(.SD, function(x) {.taperA(x)}))]
    # WoI <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
    Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)*.taperA(x)}))]
    # browser()
    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]

    names(AT) <- OCID
    names(VT) <- OCID
    names(DT) <- OCID
  }

  ## Restore Scale ----
  # browser()
  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  names(AT) <- OCID
  names(VT) <- OCID
  names(DT) <- OCID

  ## Pack Time Series - Wide Format ----
  NP <-  nrow(AT)
  ts <- seq(0,dt*(NP-1),dt)
  TSW <- data.table(ts=ts, AT = AT, VT = VT, DT = DT)
  ## Pack Time Series - Long Format ----

  IVARS <- c("ts")
  MVARS <- colnames(TSW[, -c("ts")])
  AUX <- data.table::melt(TSW, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
  TSL <- AUX[,.(t=ts,s=value,ID=gsub("\\..*$", "", variable), OCID=gsub("^[^.]*\\.", "", variable))]
  ## Trim Records, keep time series ----
  if(TrimZeros){
  TSL <- TSL[,.trimI(.SD),by=c("OCID","ID")]}

  # Rebuild AT   ----
  # browser()
  if(!is.null(Modes)){
    OCID <- colnames(AT)
    VT <- DT[, lapply(.SD, function(x) {
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    names(VT) <- OCID
    if (DetrendVT) {
      VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    }
    COLS <- colnames(AT)
    AT <- VT[, lapply(.SD, function(x) {
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    names(AT) <- OCID
  }
  if (DetrendAT) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }

  ## Return
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  return(list(TSW = TSW, TSL=TSL,Wo=Wo,PGAo=PGAo,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = UN))
}
