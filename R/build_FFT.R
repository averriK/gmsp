#' Build Fast Fourier Transform (FFT) Spectrum
#'
#' @description
#' Computes the amplitude spectrum of a time series signal using Fast Fourier
#' Transform with optional zero-padding. The function returns both single-sided
#' and double-sided amplitude spectra, with automatic frequency range limiting.
#'
#' @param .x Optional data.table containing time series data with columns 's'
#'   (signal values) and 't' (time values).
#' @param s Optional numeric vector of signal values. Required if `.x` is NULL.
#' @param t Optional numeric vector of time values. Required if `.x` is NULL.
#' @param zp Integer specifying the number of zeros to pad at the beginning and
#'   end of the signal (default: 256).
#' @param kf Numeric factor to limit the maximum frequency to kf × Nyquist
#'   frequency (default: 1.0).
#'
#' @return A data.table containing:
#'   \item{fs}{Frequency vector in Hz}
#'   \item{A}{Single-sided amplitude spectrum}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Removes the mean from the input signal
#'   \item Applies zero-padding to improve frequency resolution
#'   \item Computes FFT with length adjusted to the next power of 2
#'   \item Calculates single-sided amplitude spectrum
#'   \item Limits the frequency range based on the kf parameter
#' }
#'
#' Zero-padding improves frequency resolution but does not add new information
#' to the signal. The kf parameter allows limiting the output to frequencies
#' below the Nyquist frequency for noise reduction.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#'
#' # Create sample signal
#' dt <- 0.01 # 100 Hz sampling
#' t <- seq(0, 1, by = dt)
#' signal <- sin(2 * pi * 10 * t) + 0.5 * sin(2 * pi * 25 * t)
#'
#' # Compute FFT spectrum
#' spectrum <- build_FFT(
#'     s = signal,
#'     t = t,
#'     zp = 256,
#'     kf = 0.8
#' )
#'
#' # Plot the spectrum
#' plot(spectrum$fs, spectrum$A,
#'     type = "l",
#'     xlab = "Frequency (Hz)", ylab = "Amplitude",
#'     main = "FFT Spectrum"
#' )
#'
#' # Using data.table input
#' signal_data <- data.table(s = signal, t = t)
#' spectrum_dt <- build_FFT(.x = signal_data)
#' }
#'
#' @export
#'
#' @import data.table

build_FFT <- function(.x = NULL, s = NULL, t = NULL, zp = 256, kf = 1) {
    if (!is.null(.x)) {
        stopifnot(all(c("s", "t") %in% names(.x)))
        s <- .x$s - mean(.x$s)
        t <- .x$t
    } else {
        stopifnot(!is.null(s) & !is.null(t))
    }
    dt <- mean(diff(t))
    Fs <- 1 / dt # Sampling frequency

    # Zero-padding
    X <- c(rep(0, zp), s, rep(0, zp))
    N <- length(X)
    NP <- 2^ceiling(log2(N)) # Next power of 2 greater than N
    # Perform Fourier Transform on the zero-padded signal
    FFT <- 1 / NP * fft(X)
    # Calculate the double-sided amplitude spectra
    Ao <- Mod(FFT)
    # Calculate the single-sided amplitude spectra
    if (NP %% 2 == 0) { # NP is even
        A <- Ao[1:(NP / 2 + 1)]
        A[2:(NP / 2)] <- 2 * A[2:(NP / 2)]
        fs <- (0:(NP / 2)) * (Fs / NP)
    } else { # NP is odd
        A <- Ao[1:((NP + 1) / 2)]
        A[2:((NP - 1) / 2)] <- 2 * A[2:((NP - 1) / 2)]
        fs <- (0:((NP - 1) / 2)) * (Fs / NP)
    }
    DT <- data.table(fs = fs, A = A)

    # Cut freq
    Fnyq <- Fs / 2
    Fmax <- round(kf * Fnyq)
    DT <- DT[fs <= Fmax]
    return(DT)
}
