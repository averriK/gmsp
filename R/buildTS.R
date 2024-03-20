#' Title
#'
#' @param a data.table Time Series
#' @param dt numeric Time Step
#' @param UN character Units
#' @param Fmax integer Maximum usable frequency
#' @param FlatZeros boolean Flat Zeros Acceleration Time Series
#' @param TrimZeros boolean Trim Zeros Acceleration Time Series
#' @param PadZeros boolean Add Zeros to Acceleration Time Series
#' @param Resample boolean Add Zeros to Acceleration Time Series
#' @param RebuildAT boolean Derivate Displacements and Velocity Time Series
#' @param Detrend boolean Detrend Acceleration Time Series
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer Overlap
#'
#' @return data.table
#' @export buildTS
#'
#' @examples
#'
#' @importFrom data.table :=
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom seewave ffilter
#' @importFrom signal resample
#' @importFrom stringr str_split
#' @importFrom purrr map
#'
#'
buildTS <- function(
    a, dt, UN,
    Fmax = 25,
    Resample = TRUE,
    FlatZeros = TRUE,
    RebuildAT = FALSE,
    Detrend = TRUE,
    TrimZeros = TRUE,
    PadZeros=TRUE,
    TargetUnits = "mm",
    NW = 2048,
    OVLP = 75) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)

  . <- NULL

  OK <- is.data.table(a) && !is.null(dt) && !is.null(UN)
  stopifnot(OK)

  ATo <- copy(a)
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
  if (Detrend) {
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
  if (Detrend) {
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


  ## Integrate VT ----
  DT <- VT[, lapply(.SD, function(x) {
    x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
  })]
  names(DT) <- OCID




  # Derivate DT   ----
  if (RebuildAT) {
    COLS <- colnames(AT)
    VT <- DT[, lapply(.SD, function(x) {
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    names(VT) <- OCID

    COLS <- colnames(AT)
    AT <- VT[, lapply(.SD, function(x) {
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      return(x)
    })]
    names(AT) <- OCID
  }
  if (Detrend) {
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    # names(AT) <- OCID
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

  ## Trim Zeros
  if(TrimZeros){
    # browser()
    I <- Wo[, lapply(seq_along(.SD), function(i) {
      idx <- which(.SD[[i]]!=0) |> first()
      idx <- max(idx,2)
      return(idx)
    })] |> min()
    J <- Wo[, lapply(seq_along(.SD), function(i) {
      n <- length(.SD[[i]])
      idx <- which(.SD[[i]]!=0) |> last()
      idx <- min(idx,n-1)
      return(idx)
    })] |> max()
    AT <- AT[(I-1):(J+1)]
    VT <- VT[(I-1):(J+1)]
    DT <- DT[(I-1):(J+1)]
    Wo <- Wo[(I-1):(J+1)]
  }

  ## Restore Scale ----
  # browser()
  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * PGAo[i]})]
  names(AT) <- OCID
  names(VT) <- OCID
  names(DT) <- OCID

  ## Pack Time Series
  TS <- data.table( AT = AT, VT = VT, DT = DT)
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  NP <-  nrow(TS)
  ts <- seq(0,dt*(NP-1),dt)

  ## Return
  return(list(TS = TS, ts=ts,Wo=Wo,PGAo=PGAo,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = UN))
}
