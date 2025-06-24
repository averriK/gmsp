#' @param x          data.table.  Columns = channels, rows = samples.
#' @param ts         Numeric vector of time values.  Optional if dt is given.
#' @param dt         Numeric scalar sampling interval (seconds).  Optional if ts is given.
#' @param Units      Character.  Allowed for build_TS: "mm","cm","m","gal","g".
#'                   Allowed for build_TS_VT: "mm","cm","m".
#' @param TargetUnits Distance units after conversion ("mm","cm","m").
#' @param Fmax       Integer.  Physical maximum frequency of interest (Hz).
#' @param kNyq       Numeric (> 2.5).  Multiplier used when computing the target
#'                   resampling rate.
#' @param Resample   Logical.  If TRUE, resamples the acceleration record to
#'                   max(2.5, kNyq) * Fmax Hz.  Ignored for the raw VT input.
#' @param DetrendAT  Logical.  Subtract mean from acceleration series after each
#'                   derivation step.
#' @param DetrendVT  Logical.  Subtract mean from velocity series after each
#'                   integration or differentiation step.
#' @param DetrendDT  Logical.  Subtract mean from displacement series.
#' @param NW         Integer.  Window length for the STFT (power of two).
#' @param OVLP       Integer (0-100).  Window overlap percentage for the STFT.
#' @param FlatZeros  Logical.  If TRUE, multiplies the acceleration series by a
#'                   smooth amplitude taper.
#' @param AstopAT    Numeric.  Stop amplitude threshold for the taper (acceleration).
#' @param ApassAT    Numeric.  Pass amplitude threshold for the taper (acceleration).
#' @param TrimZeros  Logical.  If TRUE, removes rows where the taper is zero for
#'                   every channel.
#' @param Normalize  Logical.  If TRUE, divides by peak absolute amplitude before
#'                   processing and restores it afterwards.
#' @param Output     Character or NULL.  If not NULL, returns only the requested
#'                   component: "ATo","AT","VT","DT","TSW","TSL" (plus "VTo" for
#'                   build_TS_VT).  NULL returns the full list.
#' @importFrom seewave ffilter
#' @noRd
#' @export

build_TS_VT <- function(
        x, ts = NULL, dt = NULL, Units,
        DetrendVT = FALSE, DetrendAT = FALSE, DetrendDT = FALSE,
        Fmax = 16, kNyq = 3.125, Resample = TRUE,
        TargetUnits = "mm",
        NW = 128, OVLP = 75,
        FlatZeros = FALSE, AstopAT = 1e-4, ApassAT = 1e-3,
        TrimZeros = FALSE, Normalize = FALSE,
        Output = NULL) {
    
    on.exit(expr = { rm(list = ls()) }, add = TRUE)
    . <- NULL                               # for data.table
    X <- copy(x)[]
    if (!data.table::is.data.table(X)) X <- as.data.table(X)
    
    ## -------- time base ----------------------------------------------------
    if (!is.null(ts)) dt <- mean(diff(ts))
    if (is.null(dt) && is.null(ts))
        stop("Either 'dt' or 'ts' must be supplied")
    if (!is.null(dt) && !.isRational(dt)) {
        warning("Non rational dt rounded to 0.003 sec")
        dt <- round(dt, 3)
    }
    
    ## -------- unit conversion (VT) ----------------------------------------
    if (grepl(Units, "[///+]"))
        Units <- str_split(Units, "[///+]")[[1]][1]
    if (!(tolower(Units) %in% c("mm", "cm", "m")))
        stop("Units must be mm / cm / m")
    
    SFU <- if (tolower(Units) != TargetUnits)
        .getSF(SourceUnits = tolower(Units), TargetUnits = TargetUnits) else 1
    X <- X[, lapply(.SD, function(v) v * SFU)]
    OCID <- names(X)
    
    ## -------- optional normalise ------------------------------------------
    if (Normalize) {
        Vo <- apply(X, 2, function(v) max(abs(v)))
        X  <- X[, lapply(seq_along(.SD), function(i) .SD[[i]] / Vo[i])]
    }
    
    ## -------- optional detrend on VT --------------------------------------
    if (DetrendVT)
        X <- X[, lapply(.SD, function(v) v - mean(v))]
    
    ## -------- differentiate: VT -> AT -------------------------------------
    AT <- X[, lapply(.SD, .derivate, dt = dt)]
    setnames(AT, OCID)                      # keep original channel names
    
    ## -------- run the original acceleration builder -----------------------
    BT <- build_TS(
        x = AT,
        dt = dt,
        Units = TargetUnits,
        TargetUnits = TargetUnits,
        DetrendAT = DetrendAT,
        DetrendVT = DetrendVT,
        DetrendDT = DetrendDT,
        Fmax = Fmax,
        kNyq = kNyq,
        Resample = Resample,
        NW = NW, OVLP = OVLP,
        FlatZeros = FlatZeros,
        AstopAT = AstopAT, ApassAT = ApassAT,
        TrimZeros = TrimZeros,
        Normalize = FALSE,              # done here if requested
        Output = NULL
    )
    
    ## -------- restore magnitude if normalised -----------------------------
    if (Normalize) {
        mult <- function(tab, vec) {
            tab[, lapply(seq_along(.SD), function(i) .SD[[i]] * vec[i])]
        }
        BT$AT <- mult(BT$AT, Vo)
        BT$VT <- mult(BT$VT, Vo)
        BT$DT <- mult(BT$DT, Vo)
        
        # wide table
        for (i in seq_along(OCID)) {
            ch <- OCID[i]; v <- Vo[i]
            BT$TSW[[paste0("AT.", ch)]] <- BT$TSW[[paste0("AT.", ch)]] * v
            BT$TSW[[paste0("VT.", ch)]] <- BT$TSW[[paste0("VT.", ch)]] * v
            BT$TSW[[paste0("DT.", ch)]] <- BT$TSW[[paste0("DT.", ch)]] * v
        }
        
        # long table
        scaleDT <- data.table(OCID = OCID, v = Vo)
        BT$TSL <- merge(BT$TSL, scaleDT, by = "OCID")
        BT$TSL[, s := s * v][, v := NULL]
    }
    
    ## -------- fix metadata -------------------------------------------------
    BT$SourceUnits <- Units              # keep original velocity units
    # list names and structure are otherwise unchanged
    
    ## -------- selector identical to build_TS ------------------------------
    if (!is.null(Output)) {
        if (!(Output %in% c("ATo", "AT", "VT", "DT", "TSW", "TSL")))
            stop("Unknown Output flag")
        return(BT[[Output]])
    }
    
    return(BT)
}



