#' Extract Intrinsic Mode Functions (IMFs) from Time Series
#'
#' @description
#' Decomposes a time series signal into Intrinsic Mode Functions (IMFs) using
#' Empirical Mode Decomposition (EMD), Ensemble EMD (EEMD), or Complete Ensemble
#' EMD (CEEMD). IMFs represent oscillatory modes embedded in the signal, each
#' with a characteristic frequency and amplitude.
#'
#' @param .x Optional data.table containing time series data with columns 's'
#'   (signal values) and 't' (time values).
#' @param s Optional numeric vector of signal values. Required if `.x` is NULL.
#' @param ts Optional numeric vector of time values. Required if `.x` is NULL.
#' @param dt Optional numeric value for time step. If NULL, calculated from `ts` vector.
#' @param method Character string specifying the decomposition method:
#'   \itemize{
#'     \item "emd" - Empirical Mode Decomposition (default)
#'     \item "eemd" - Ensemble Empirical Mode Decomposition
#'     \item "ceemd" - Complete Ensemble Empirical Mode Decomposition
#'   }
#' @param boundary Character string specifying the boundary condition:
#'   \itemize{
#'     \item "wave" - Wave boundary condition (default)
#'     \item "none" - No boundary condition
#'     \item "symmetric" - Symmetric boundary condition
#'     \item "periodic" - Periodic boundary condition
#'     \item "evenodd" - Even-odd boundary condition
#'   }
#' @param max.imf Integer specifying the maximum number of IMFs to extract (default: 15).
#' @param noise.type Character string for noise type in ensemble methods:
#'   \itemize{
#'     \item "gaussian" - Gaussian noise (default)
#'     \item "uniform" - Uniform noise
#'   }
#' @param noise.amp Numeric value for noise amplitude in ensemble methods (default: 0.5e-7).
#' @param trials Integer specifying the number of trials for ensemble methods (default: 10).
#' @param stop.rule Character string for the stopping rule:
#'   \itemize{
#'     \item "type5" - Type 5 stopping rule (default)
#'     \item "type1" - Type 1 stopping rule
#'     \item "type2" - Type 2 stopping rule
#'     \item "type3" - Type 3 stopping rule
#'     \item "type4" - Type 4 stopping rule
#'   }
#' @param verbose Logical. If TRUE, prints progress information for ensemble methods.
#'
#' @return A data.table in long format containing:
#'   \item{t}{Time values}
#'   \item{s}{Signal values (IMFs, original signal, and residue)}
#'   \item{IMF}{Component identifier (IMF1, IMF2, ..., signal, residue)}
#'
#' @details
#' The function implements three different EMD variants:
#' \enumerate{
#'   \item \strong{EMD}: Standard Empirical Mode Decomposition that iteratively
#'         extracts IMFs by sifting the signal
#'   \item \strong{EEMD}: Ensemble EMD adds white noise to the signal multiple
#'         times and averages the results to reduce mode mixing
#'   \item \strong{CEEMD}: Complete Ensemble EMD uses complementary noise pairs
#'         to further improve the decomposition quality
#' }
#'
#' Each IMF represents a mode of oscillation with varying amplitude and frequency,
#' making it useful for analyzing non-stationary and non-linear signals.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#'
#' # Create sample signal with multiple frequencies
#' dt <- 0.01 # 100 Hz sampling
#' t <- seq(0, 2, by = dt)
#' signal <- sin(2 * pi * 2 * t) + 0.5 * sin(2 * pi * 10 * t) + 0.1 * rnorm(length(t))
#'
#' # Extract IMFs using standard EMD
#' imfs <- get_imf(
#'     s = signal,
#'     ts = t,
#'     method = "emd",
#'     max.imf = 5
#' )
#'
#' # Extract IMFs using ensemble EMD
#' imfs_eemd <- get_imf(
#'     s = signal,
#'     ts = t,
#'     method = "eemd",
#'     trials = 20,
#'     noise.amp = 0.1
#' )
#'
#' # Plot the first few IMFs
#' library(ggplot2)
#' ggplot(
#'     imfs[IMF %in% c("IMF1", "IMF2", "signal")],
#'     aes(x = t, y = s, color = IMF)
#' ) +
#'     geom_line() +
#'     facet_wrap(~IMF, ncol = 1, scales = "free_y")
#' }
#'
#' @export
#'
#' @importFrom data.table :=
#' @importFrom stats na.omit
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom hht EEMDCompile
#' @importFrom hht CEEMD
#' @importFrom hht EEMD
#' @importFrom EMD emd
#'
get_imf <- function(.x = NULL, s = NULL, ts = NULL, dt = NULL, method = "emd", boundary = "wave", max.imf = 15, noise.type = "gaussian", noise.amp = 0.5e-7, trials = 10, stop.rule = "type5", verbose = FALSE) {
    on.exit(expr = {
        rm(list = ls())
    }, add = TRUE)

    stopifnot(tolower(method) %in% c("emd", "eemd", "ceemd"))
    stopifnot(tolower(stop.rule) %in% c("type1", "type2", "type3", "type4", "type5"))
    stopifnot(tolower(boundary) %in% c("none", "wave", "symmetric", "periodic", "evenodd"))
    stopifnot(tolower(noise.type) %in% c("uniform", "gaussian"))

    if (is.null(s) && is.null(.x)) {
        stop("No acceleration data provided.")
    }
    if (!is.null(.x)) {
        s <- .x$s
        ts <- .x$t
    }
    if (is.null(ts) && is.null(.x) && is.null(dt)) {
        stop("No time data provided.")
    }

    if (!is.null(ts)) {
        dt <- mean(diff(ts))
        to <- min(ts)
    }

    ts <- seq(0, (length(s) - 1) * dt, by = dt)

    DT <- data.table(t = ts, s = s)
    DT <- .trimZeros(DT)
    n <- nrow(DT)
    rm(s)

    if (tolower(method) == "emd") {
        AUX <- EMD::emd(xt = DT$s, tt = DT$t, boundary = boundary, max.imf = max.imf, stoprule = stop.rule)
        IMF <- as.data.table(AUX$imf)
        RES <- AUX$residue
    }

    if (tolower(method) == "eemd") {
        nimf <- EMD::emd(xt = DT$s, tt = DT$t, boundary = boundary, max.imf = max.imf, stoprule = stop.rule)$nimf
        DIR <- tempdir(check = TRUE)
        hht::EEMD(sig = DT$s, tt = DT$t, nimf = nimf, max.imf = max.imf, boundary = boundary, noise.amp = noise.amp, noise.type = noise.type, trials = trials, stop.rule = stop.rule, trials.dir = DIR, verbose = verbose)
        AUX <- EEMDCompile(trials.dir = DIR, trials = trials, nimf = nimf) |> suppressWarnings()
        unlink(DIR, force = TRUE, recursive = TRUE)
        IMF <- as.data.table(AUX$averaged.imfs)
        RES <- AUX$averaged.residue
    }

    if (tolower(method) == "ceemd") {
        AUX <- hht::CEEMD(sig = DT$s, tt = DT$t, noise.amp = noise.amp, noise.type = noise.type, trials = trials, stop.rule = stop.rule, verbose = verbose)
        IMF <- as.data.table(AUX$imf)
        RES <- AUX$residue
    }

    if (ncol(IMF) == 0) {
        warning("No IMFs were found.")
        return(NULL)
    }
    # Zero padding to ensure IMF has the same length as s and t
    if (nrow(IMF) < nrow(DT)) {
        pad_length <- nrow(DT) - nrow(IMF)
        padding <- data.table(matrix(0, nrow = pad_length, ncol = ncol(IMF)))
        setnames(padding, names(IMF))
        IMF <- rbind(IMF, padding)
    }

    # Rename columns before adding signal and residue
    setnames(IMF, c(paste0("IMF", seq_len(ncol(IMF)))))

    # Add signal and residue as the last columns of IMF
    TSW <- data.table(t = DT$t + to, signal = DT$s, IMF, residue = RES)
    #
    ivars <- c("t")
    mvars <- colnames(TSW[, -c("t")])
    AUX <- melt(TSW, id.vars = ivars, measure.vars = mvars)
    TSL <- AUX[, list(t = t, s = value, IMF = variable)]

    return(TSL)
}
