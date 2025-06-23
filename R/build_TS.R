# nolint start
#' Build Time Series from Seismic Data
#'
#' @description
#' Processes seismic acceleration time series data to generate acceleration (AT),
#' velocity (VT), and displacement (DT) time series. The function performs unit
#' conversion, resampling, filtering, integration/differentiation, and optional
#' detrending and normalization.
#'
#' @param x A data.table containing seismic acceleration time series data.
#'   Each column represents a different component or station.
#' @param ts Optional numeric vector of time values. If provided, `dt` will be
#'   calculated as the mean time step.
#' @param dt Optional numeric value for the time step (sampling interval) in seconds.
#' @param Units Character string specifying the units of the input data.
#'   Must be one of: "mm", "cm", "m", "gal", "g". If the string contains
#'   separators like "///" or "+", only the first part is used.
#' @param DetrendAT Logical. If TRUE, removes the mean from acceleration time series.
#' @param DetrendVT Logical. If TRUE, removes the mean from velocity time series.
#' @param DetrendDT Logical. If TRUE, removes the mean from displacement time series.
#' @param Fmax Integer specifying the maximum frequency for filtering (default: 16 Hz).
#' @param kNyq Numeric value for Nyquist frequency factor (default: 3.125, must be > 2.5).
#' @param Resample Logical. If TRUE, resamples the data to a target frequency
#'   calculated as max(2.5, kNyq) * Fmax (default: TRUE).
#' @param TargetUnits Character string for the target units after conversion
#'   (default: "mm").
#' @param NW Integer specifying the window length for FFT operations (default: 128).
#'   Will be adjusted to the nearest power of 2 if it exceeds half the data length.
#' @param OVLP Integer specifying the overlap percentage for FFT windows (default: 75).
#' @param FlatZeros Logical. If TRUE, applies amplitude-based tapering to flatten
#'   near-zero regions (default: FALSE).
#' @param AstopAT Numeric value for the stop amplitude threshold in acceleration
#'   tapering (default: 1e-4).
#' @param ApassAT Numeric value for the pass amplitude threshold in acceleration
#'   tapering (default: 1e-3).
#' @param TrimZeros Logical. If TRUE, removes time points where all components
#'   have zero amplitude (default: FALSE).
#' @param Normalize Logical. If TRUE, normalizes the data by dividing by the
#'   maximum absolute value of each component (default: FALSE).
#' @param Output Character string specifying the output format. If NULL (default),
#'   returns a complete list. Options: "ATo", "AT", "VT", "DT", "TSW", "TSL".
#'
#' @return A list containing the processed time series data and metadata:
#'   \item{ATo}{Original acceleration time series scaled to target units}
#'   \item{TSL}{Long format time series data with columns: t (time), s (value),
#'              ID (component type), OCID (original column name)}
#'   \item{TSW}{Wide format time series data with columns: ts (time), AT, VT, DT}
#'   \item{Wo}{Tapering window applied to the data}
#'   \item{Fs}{Final sampling frequency in Hz}
#'   \item{dt}{Final time step in seconds}
#'   \item{df}{Frequency resolution in Hz}
#'   \item{fs}{Frequency vector for spectral analysis}
#'   \item{NP}{Number of data points}
#'   \item{TargetUnits}{Target units used for conversion}
#'   \item{SourceUnits}{Original units of input data}
#'
#'   If `Output` is specified, returns only the requested component:
#'   \item{ATo}{Original acceleration time series}
#'   \item{AT}{Processed acceleration time series}
#'   \item{VT}{Velocity time series}
#'   \item{DT}{Displacement time series}
#'   \item{TSW}{Wide format time series}
#'   \item{TSL}{Long format time series}
#'
#' @details
#' The function performs the following processing steps:
#' \enumerate{
#'   \item Unit conversion to target units if necessary
#'   \item Optional normalization by maximum absolute value
#'   \item Optional detrending by removing mean
#'   \item Resampling to target frequency if enabled
#'   \item Integration of acceleration to velocity using frequency-domain filtering
#'   \item Integration of velocity to displacement
#'   \item Differentiation back to velocity and acceleration for consistency
#'   \item Optional amplitude-based tapering for near-zero regions
#'   \item Optional trimming of zero-amplitude regions
#' }'
#' @export
#'
#' @import data.table
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom seewave ffilter
#' @importFrom signal resample
#' @importFrom stringr str_split
#' @importFrom purrr map
#'
build_TS <- function(
    x, ts = NULL, dt = NULL, Units,
    DetrendAT = FALSE,
    DetrendVT = FALSE,
    DetrendDT = FALSE,
    Fmax = 16,
    kNyq = 3.125, #>2.5
    Resample = TRUE,
    TargetUnits = "mm",
    NW = 128,
    OVLP = 75,
    FlatZeros = FALSE,
    AstopAT = 1e-4,
    ApassAT = 1e-3,
    TrimZeros = FALSE,
    Normalize = FALSE,
    Output = NULL) {
    X <- copy(x) |> as.data.table()

    NP <- nrow(X)
    NW_min <- NP / 2
    if (NW > NW_min) {
        L2 <- 2^floor(log2(NW_min))
        U2 <- 2^ceiling(log2(NW_min))
        if (abs(L2 - NW_min) < abs(U2 - NW_min)) {
            NW <- min(as.integer(L2), NP)
        } else {
            NW <- min(as.integer(U2), NP)
        }
    }
    if (!is.null(ts)) {
        dts <- diff(ts)
        dt <- mean(dts)
    }

    # Check if dt is a rational number
    if (!.isRational(dt)) {
        warning("Time step is not a rational number. Rounding to 3 decimal places (max sampling frequency 1kHz)")
        dt <- round(dt, 3)
    }
    ts <- seq(0, (NP - 1) * dt, by = dt)
    Fs <- 1 / dt

    ## Scale Units  ----
    SFU <- 1

    if (grepl(Units, pattern = "[///+]")) {
        Units <- (str_split(Units, pattern = "[///+]") |> unlist())[1]
    }
    if (!(tolower(Units) %in% c("mm", "cm", "m", "gal", "g"))) {
        return(NULL)
    }
    if (tolower(Units) != TargetUnits) {
        SFU <- .getSF(SourceUnits = tolower(Units), TargetUnits = TargetUnits)
        X <- X[, .(sapply(.SD, function(x) {
            x * SFU
        }))]
    }

    OCID <- names(X)

    # Export RAW record scaled to TargetUnits

    ATo <- data.table(ts = ts, Units = TargetUnits, X)
    if (!is.null(Output) && Output == "ATo") {
        return(ATo)
    }

    # setnames(RTSW,old=OCID,new=paste0("AT.",OCID))

    ## Scale record ----
    if (Normalize == TRUE) {
        Ao <- apply(X, 2, function(x) {
            max(abs(x))
        })
        X <- X[, .(sapply(.SD, function(x) {
            x / max(abs(x))
        }))]
    }

    # Detrend
    if (DetrendAT == TRUE) {
        X <- X[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }


    ## Resample ----
    if (Resample) {
        TargetFs <- as.integer(max(2.5, kNyq) * Fmax) # 80/200 Hz
        X <- .resample(X, dt = dt, TargetFs = TargetFs, Fmax = Fmax, NW = NW, OVLP = OVLP)
        Fs <- TargetFs
        dt <- 1 / Fs
    }
    ## Case #2. Acceleration Time Histories ----
    if (FlatZeros == TRUE) {
        Wo <- X[, .(sapply(.SD, function(x) {
            .taperA(x, Astop = SFU * AstopAT, Apass = SFU * ApassAT)
        }))]
    } else {
        Wo <- X[, .(sapply(.SD, function(x) {
            rep(1, length(x))
        }))]
    }

    AT <- X
    AT <- AT[, lapply(seq_along(.SD), function(i) {
        .SD[[i]] * Wo[[i]]
    })]
    # AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    # Detrend
    if (DetrendAT == TRUE) {
        AT <- AT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }


    names(AT) <- OCID
    # browser()
    VT <- AT[, lapply(.SD, function(x) {
        .integrate(dx = x, dt = dt, NW = NW, OVLP = OVLP)
    })]

    if (nrow(VT) > nrow(Wo)) {
        O <- data.table(sapply(Wo, function(x) rep(0, nrow(VT) - nrow(Wo))))
        Wo <- rbind(Wo, O)
    }
    VT <- VT[, lapply(seq_along(.SD), function(i) {
        .SD[[i]] * Wo[[i]]
    })]
    # VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    if (DetrendVT == TRUE) {
        VT <- VT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }


    names(VT) <- OCID

    DT <- VT[, lapply(.SD, function(x) {
        .integrate(dx = x, dt = dt, NW = NW, OVLP = OVLP)
    })]

    if (nrow(DT) > nrow(Wo)) {
        O <- data.table(sapply(Wo, function(x) rep(0, nrow(DT) - nrow(Wo))))
        Wo <- rbind(Wo, O)
    }
    DT <- DT[, lapply(seq_along(.SD), function(i) {
        .SD[[i]] * Wo[[i]]
    })]

    ## Detrend DT ----
    if (DetrendDT == TRUE) {
        DT <- DT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }
    names(DT) <- OCID




    ## Derivate DT ----
    VT <- DT[, lapply(.SD, function(x) {
        .derivate(X = x, dt = dt)
    })]

    if (DetrendVT == TRUE) {
        VT <- VT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }
    ## Derivate VT ----
    AT <- VT[, lapply(.SD, function(x) {
        .derivate(X = x, dt = dt)
    })]

    ## Detrend AT ----
    if (DetrendAT == TRUE) {
        AT <- AT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }
    ## Taper Zeros ----
    # Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]

    if (FlatZeros == TRUE) {
        Wo <- AT[, .(sapply(.SD, function(x) {
            .taperA(x, Astop = SFU * AstopAT, Apass = SFU * ApassAT)
        }))]
    } else {
        Wo <- AT[, .(sapply(.SD, function(x) {
            rep(1, length(x))
        }))]
    }


    AT <- AT[, lapply(seq_along(.SD), function(i) {
        .SD[[i]] * Wo[[i]]
    })]
    VT <- VT[, lapply(seq_along(.SD), function(i) {
        .SD[[i]] * Wo[[i]]
    })]
    DT <- DT[, lapply(seq_along(.SD), function(i) {
        .SD[[i]] * Wo[[i]]
    })]


    # Trim Zeros
    if (TrimZeros) {
        idx <- apply(Wo != 0, MARGIN = 1, function(x) {
            all(x)
        })
        AT <- AT[idx]
        VT <- VT[idx]
        DT <- DT[idx]
    }

    # Fix trend
    if (DetrendAT == TRUE) {
        AT <- AT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }
    if (DetrendVT == TRUE) {
        VT <- VT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }
    if (DetrendDT == TRUE) {
        DT <- DT[, .(sapply(.SD, function(x) {
            x - mean(x)
        }))]
    }



    ## Normalize Maximum ----
    if (Normalize == TRUE) {
        AT <- AT[, lapply(seq_along(.SD), function(i) {
            .SD[[i]] * Ao[i]
        })]
        VT <- VT[, lapply(seq_along(.SD), function(i) {
            .SD[[i]] * Ao[i]
        })]
        DT <- DT[, lapply(seq_along(.SD), function(i) {
            .SD[[i]] * Ao[i]
        })]
    }

    names(AT) <- OCID
    names(VT) <- OCID
    names(DT) <- OCID
    if (!is.null(Output) && Output == "AT") {
        return(AT)
    }
    if (!is.null(Output) && Output == "VT") {
        return(VT)
    }
    if (!is.null(Output) && Output == "DT") {
        return(DT)
    }

    ## Summary ----
    NP <- nrow(AT)
    ts <- seq(0, dt * (NP - 1), dt)
    Fs <- 1 / dt
    df <- Fs / NW # 0.03125#
    fs <- seq(from = 0, by = df, length.out = NW / 2)
    ## Pack Time Series  ----

    TSW <- data.table(ts = ts, AT = AT, VT = VT, DT = DT)
    if (!is.null(Output) && Output == "TSW") {
        return(TSW)
    }

    ivars <- c("ts")
    mvars <- colnames(TSW[, -c("ts")])
    AUX <- data.table::melt(TSW, id.vars = ivars, measure.vars = mvars) |> na.omit()
    TSL <- AUX[, .(t = ts, s = value, ID = gsub("\\..*$", "", variable), OCID = gsub("^[^.]*\\.", "", variable))]
    if (!is.null(Output) && Output == "TSL") {
        return(TSL)
    }

    ## Trim Records,  ----
    # if(TrimZeros){
    #   #Warning: Different time scales from now on
    #   TSL <- TSL[,.trimZeros(.SD),by=c("OCID","ID")]}

    ## Return ----
    return(list(ATo = ATo, TSL = TSL, TSW = TSW, Wo = Wo, Fs = Fs, dt = dt, df = df, fs = fs, NP = NP, TargetUnits = TargetUnits, SourceUnits = Units))
}
