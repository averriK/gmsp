#' Title
#'
#' @param a data.table Time Series
#' @param dt numeric Time Step
#' @param UN character Units
#' @param hl numeric 1/Tolerance
#' @param DownFs integer Downsample Frequency
#' @param UpFs integer Upsample Frequency
#' @param FlatZerosAT boolean Flat Zeros Acceleration Time Series
#' @param DerivateDT boolean Derivate Displacements Time Series
#' @param DerivateVT boolean Derivate Velocity Time Series
#' @param DetrendAT boolean Detrend Acceleration Time Series
#' @param DetrendVT boolean Detrend Velocity Time Series
#' @param DetrendDT boolean Detrend Displacement Time Series
#' @param Fpass_LP numeric Low Pass Frequency
#' @param Fstop_LP numeric Low Stop Frequency
#' @param Fpass_HP numeric High Pass Frequency
#' @param Fstop_HP numeric High Stop Frequency
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer Overlap
#'
#' @return data.table
#' @export buildTS
#'
#' @examples
#'
#' @import data.table
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom seewave ffilter
#' @importFrom signal resample
#' @importFrom pracma detrend
#' @importFrom stringr str_split
#' @importFrom purrr map
#'
#'
buildTS <- function(
    a, dt, UN = NULL, hl = 500, # AT<PGAo/500->0
    DownFs = 0, UpFs = 0, FlatZerosAT = TRUE, DerivateDT = FALSE, DerivateVT = FALSE,
    DetrendAT = TRUE, DetrendVT = FALSE, DetrendDT = FALSE,
    Fpass_LP = 0, Fstop_LP = 0, Fpass_HP = 0, Fstop_HP = 0,
    TargetUnits = "mm", NW = 2048, OVLP = 75) {
  on.exit(expr = {
    rm(list = ls())
  }, add = TRUE)


  OK <- is.data.table(a) # && !is.null(dt) && !is.null(UN)
  stopifnot(OK)

  ATo <- copy(a)
  ## Check Record ----
  NP <- nrow(ATo)
  if (  NP == 0 ||dt == 0 ||any(is.na(ATo)) ||max(abs(ATo)) == 0 ) {
    # Null Record
    return(NULL)
  }
  ## Check Length ----
  if (NP < NW) return(NULL)
  OCID <- names(ATo)
  ## Check Invalid Header ----
  if (length(unique(OCID)) < 3) return(NULL)


  ## Scale Units ----

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
  ## Scale to 1 ----
  AT <-ATo[, .(sapply(.SD, function(x){x/max(abs(x))}))]


  ## Detrend ----
  if (DetrendAT) {
    AT[, (colnames(AT)) := lapply(.SD, function(x) {
      pracma::detrend(x, tt = "linear")
    })]
  }


  ## Flat Zeros & Taper  ----
  if (FlatZerosAT == TRUE) {
    # browser()
    AT <- AT[, lapply(X = .SD, FUN = function(x) {
      n <- length(x)
      Astop <- max(PGAo / hl, min(x))
      Apass <- max(PGAo / hl / 2, min(x))

      iH_stop <- which(abs(x) > Astop) |> first() # AT<PGAo/hl->stop
      iL_stop <- which(abs(x) > Astop) |> last()
      iH_pass <- which(abs(x) > Apass) |> first() # AT<PGAo/hl/2->pass
      iL_pass <- which(abs(x) > Apass) |> last()

      if (length(iH_pass) == 1 && length(iH_stop) == 1 && iH_pass > iH_stop) {
        HP <- .buildHighPassButtterworth(f = seq(1, n), Fpass = iH_pass, Fstop = iH_stop, Astop = 0.001, Apass = 0.999)
      } else {
        HP <- rep(1, n)
      }
      if (length(iL_pass) == 1 && length(iL_stop) == 1 && iL_pass < iL_stop) {
        LP <- .buildLowPassButtterworth(f = seq(1, n), Fstop = iL_stop, Fpass = iL_pass, Astop = 0.001, Apass = 0.999)
      } else {
        LP <- rep(1, n)
      }
      return(x * LP * HP)
    })]
  }



  ## Up-sampling ----

  if (UpFs != 0 && UpFs > Fs) {
    AT <- AT[, lapply(.SD, function(x) {
      x <- signal::resample(x, UpFs, Fs)
      # x <- pracma::detrend(x,tt="linear")
      return(x)
    })]

    ## Detrend ----
    if (DetrendAT) {
      AT[, (colnames(AT)) := lapply(.SD, function(x) {
        pracma::detrend(x, tt = "linear")
      })]
    }

    # Update Record
    Fs <- UpFs
    dt <- 1 / Fs
    NP <- nrow(AT)
    names(AT) <- OCID
  }

  ## Padding Zeros ----
  NP <- nrow(AT)
  NZ <- .getNZ(NP)
  if (NZ > 0) {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NZ))]
  } else {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NW))]
  }
  AT <- rbindlist(list(ZEROS, AT))



  ## Antialias & Downsampling ----

  if (DownFs != 0 && Fs != DownFs) {
    df <- Fs / NW # 0.03125#
    fs <- seq(from = 0, by = df, length.out = NW / 2)
    FsNYQ <- Fs / 2 # 16
    COLS <- colnames(AT)

    if (Fstop_LP == 0 & Fpass_LP == 0) {
      LP <- rep(1, times = NW / 2)
    } else {
      LP <- .buildLowPassButtterworth(f = fs, Fstop = round(1 * Fstop_LP / df) * df, Fpass = round(1 * Fpass_LP / df) * df, Astop = 0.001, Apass = 0.95)
    }


    AT <- AT[, lapply(.SD, function(x) {
      x <- ffilter(wave = x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      x <- signal::resample(x, DownFs, Fs)
      return(x)
    })]

    ## Detrend ---
    if (DetrendAT) {
      AT <- AT[, lapply(.SD, function(x) {
        x <- pracma::detrend(x, tt = "linear")
        return(x)
      })]
    }
    names(AT) <- COLS
    # Update
    Fs <- DownFs
    dt <- 1 / Fs
  }

  ## Padding zeros ----

  NP <- nrow(AT)
  if (NP < NW) {
    return(NULL)
  }
  NZ <- .getNZ(NP)
  if (NZ > 0) {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NZ))]
  } else {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NW))]
  }
  AT <- rbindlist(list(ZEROS, AT))

  ## Build HT Filters ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  AT <- as.data.table(AT)






  ## Padding zeros ----
  NP <- nrow(AT)
  NZ <- .getNZ(NP)
  if (NZ > 0) {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NZ))]
  } else {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NW))]
  }
  AT <- rbindlist(list(ZEROS, AT))
  ## Build Filters ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)


  FsNYQ <- Fs / 2 # 30
  if (Fstop_LP == 0 & Fpass_LP == 0) {
    LP <- rep(1, times = NW / 2)
  } else {
    LP <- .buildLowPassButtterworth(f = fs, Fstop = round(1 * Fstop_LP / df) * df, Fpass = round(1 * Fpass_LP / df) * df, Astop = 0.001, Apass = 0.90)
  }

  #
  if (Fpass_HP == 0 & Fstop_HP == 0) {
    HP <- rep(1, times = NW / 2)
  } else {
    HP <- .buildHighPassButtterworth(f = fs, Fstop = round(1 * Fstop_HP / df) * df, Fpass = round(1 * Fpass_HP / df) * df, Astop = 0.001, Apass = 0.99)
  }
  ## Highpass & Lowpass filters----
  COLS <- colnames(AT)
  AT <- AT[, lapply(.SD, function(x) {
    x <- ffilter(wave = x, f = Fs, wl = NW, ovlp = OVLP, custom = HP * LP, rescale = TRUE)
    # x <- pracma::detrend(x,tt="linear")
    return(x)
  })]
  ## Detrend ----
  if (DetrendAT) {
    AT <- AT[, lapply(.SD, function(x) {
      x <- pracma::detrend(x, tt = "linear")
      return(x)
    })]
  }
  names(AT) <- COLS

  ## Build VT Filters ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  FsNYQ <- Fs / 2 # 16
  if (Fstop_LP == 0 & Fpass_LP == 0) {
    LP <- rep(1, times = NW / 2)
  } else {
    LP <- .buildLowPassButtterworth(f = fs, Fstop = round(Fstop_LP / df) * df, Fpass = round(Fpass_LP / df) * df, Astop = 0.001, Apass = 0.90)
  }

  HI <- .buildIntegrateFilter(f = fs) ## Integrate Filter

  ## Padding zeros AT----
  NP <- nrow(AT)
  NZ <- .getNZ(NP)
  if (NZ > 0) {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NZ))]
  } else {
    ZEROS <- data.table()[, (colnames(AT)) := list(rep(0, NW))]
  }
  AT <- rbindlist(list(ZEROS, AT))
  # Integrate Acceleration ----
  COLS <- colnames(AT)
  VT <- AT[, lapply(.SD, function(x) {
    x <- .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
    # x <- pracma::detrend(x,tt="linear")
  })]
  names(VT) <- COLS

  ## Build Displacement Filters ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  FsNYQ <- Fs / 2 # 16 round(3/8*Fpass_LP/df)*df
  if (Fstop_LP == 0 & Fpass_LP == 0) {
    LP <- rep(1, times = NW / 2)
  } else {
    LP <- .buildLowPassButtterworth(f = fs, Fstop = round(Fstop_LP / df) * df, Fpass = round(Fpass_LP / df) * df, Astop = 0.001, Apass = 0.90)
  }
  HI <- .buildIntegrateFilter(f = fs) ## Integrate Filter

  ## Padding zeros VT ----
  NP <- nrow(VT)
  NZ <- .getNZ(NP)
  if (NZ > 0) {
    ZEROS <- data.table()[, (colnames(VT)) := list(rep(0, NZ))]
  } else {
    ZEROS <- data.table()[, (colnames(VT)) := list(rep(0, NW))]
  }
  VT <- rbindlist(list(ZEROS, VT))


  ## Integrate Velocity ----
  COLS <- colnames(AT)
  DT <- VT[, lapply(.SD, function(x) {
    x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
    # x <- pracma::detrend(x,tt="linear")
  })]
  names(DT) <- COLS

  ## Build Derivative Filters ----
  FsNYQ <- Fs / 2 # 30

  if (Fstop_LP == 0 & Fpass_LP == 0) {
    LP <- rep(1, times = NW / 2)
  } else {
    LP <- .buildLowPassButtterworth(f = fs, Fstop = round(1 * Fstop_LP / df) * df, Fpass = round(1 * Fpass_LP / df) * df, Astop = 0.001, Apass = 0.90)
  }

  HD <- .buildDerivateFilter(f = fs) ## Derivate Filter

  if (Fpass_HP == 0 & Fstop_HP == 0) {
    HP <- rep(1, times = NW / 2)
  } else {
    HP <- .buildHighPassButtterworth(f = fs, Fstop = round(1 * Fstop_HP / df) * df, Fpass = round(1 * Fpass_HP / df) * df, Astop = 0.001, Apass = 0.99)
  }

  # HP <-  .buildHighPassButtterworth(f=fs, Fstop = Fstop_HP,Fpass=Fpass_HP ,Astop=0.001, Apass=0.99)



  # Derivate DT   ----
  if (DerivateDT) {
    if (Fstop_LP == 0 & Fpass_LP == 0) {
      LP <- rep(1, times = NW / 2)
    } else {
      LP <- .buildLowPassButtterworth(f = fs, Fstop = round(Fstop_LP / df) * df, Fpass = round(Fpass_LP / df) * df, Astop = 0.001, Apass = 0.90)
    }

    COLS <- colnames(AT)
    VT <- DT[, lapply(.SD, function(x) {
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      # x <- pracma::detrend(x,tt="linear")
      return(x)
    })]
    names(VT) <- COLS
  }


  # Derivate VT   ----
  if (DerivateVT) {
    COLS <- colnames(AT)
    AT <- VT[, lapply(.SD, function(x) {
      x <- NW * .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HP * LP, rescale = TRUE)
      # x <- pracma::detrend(x,tt="linear")
      return(x)
    })]
    names(AT) <- COLS
  }






  ## Pack & Taper  ----

  NMX <- min(nrow(AT), nrow(VT), nrow(DT))
  AT <- AT[-((NMX):.N)]
  VT <- VT[-((NMX):.N)]
  DT <- DT[-((NMX):.N)]
  PGA <- apply(AT, 2, function(x) {
    max(abs(x))
  })
  Wo <- .taper(AT[[which.max(PGA)]])

  AT <- AT[,.(sapply(.SD, function(x) {x * Wo}))]
  VT <- VT[,.(sapply(.SD, function(x) {x * Wo}))]
  DT <- DT[,.(sapply(.SD, function(x) {x * Wo}))]


  ## Restore Scale ----
  browser()
  ATo <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * PGAo[i])]
  VTo <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * PGAo[i])]
  DTo <- AT[, lapply(seq_along(.SD), function(i) .SD[[i]] * PGAo[i])]

  ## Pack Time Series
  TS <- list(data.table(W = Wo, AT = ATo, VT = VTo, DT = DTo))
  names(TS) <- "I"


  ## Trim Zeros ---
  TS <- lapply(TS, unique)
  ## Add row of Zeros
  TS$I <- rbindlist(list(TS$I, 0 * TS$I[1]))

  ## Return
  return(list(TS = TS, Fs = Fs, dt = dt, NP = NP, UN = "mm"))
}
