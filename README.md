# Ground Motion Signal Processing Tools (gmsp)

[![R-CMD-check](https://github.com/your-username/gmsp/workflows/R-CMD-check/badge.svg)](https://github.com/your-username/gmsp/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/gmsp)](https://CRAN.R-project.org/package=gmsp)
[![License](https://img.shields.io/badge/license-file%20LICENSE-blue.svg)](LICENSE)

## Overview

The `gmsp` package provides a comprehensive set of tools for processing ground motion signals, specifically designed for earthquake engineering and seismological applications. It includes functions for spectral analysis, time-frequency decomposition, intensity measure calculations, and advanced signal processing techniques including the Hilbert-Huang Transform.

## Features

### Core Signal Processing
- **Response Spectra Calculation**: Compute pseudo-spectral acceleration using state-space formulation
- **FFT Analysis**: Fast Fourier Transform with zero-padding for improved frequency resolution
- **Spectral Smoothing**: Multiple smoothing algorithms (moving average, Savitzky-Golay, Gaussian kernel, etc.)

### Advanced Decomposition Methods
- **Empirical Mode Decomposition (EMD)**: Extract intrinsic mode functions from non-stationary signals
- **Ensemble EMD (EEMD)**: Noise-assisted EMD for improved mode separation
- **Complete Ensemble EMD (CEEMD)**: Enhanced EEMD with complementary noise pairs

### Intensity Measures
- **Peak Ground Motion**: PGA, PGV, PGD calculations
- **Arias Intensity**: Total energy content of ground motion
- **Duration Measures**: Bracketed (Hussid) duration calculations
- **Damage Indices**: PPI, EPI, PDI for structural damage assessment

### Specialized Tools
- **Newmark Displacement**: Sliding displacement analysis for slope stability
- **Time Series Processing**: Filtering, baseline correction, and data manipulation
- **Multi-format Support**: Handle various accelerogram formats and data structures

## Installation

### From CRAN (when available)
```r
install.packages("gmsp")
```

### Development Version
```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install gmsp from GitHub
devtools::install_github("your-username/gmsp")
```

### Dependencies
The package requires R >= 4.0 and depends on several specialized packages:
```r
# Core dependencies
install.packages(c("data.table", "pracma", "seewave", "signal", "zoo", 
                   "stats", "stringr", "purrr", "digest", "spectral", 
                   "hht", "EMD", "expm"))
```

## Quick Start

### Load the package
```r
library(gmsp)
library(data.table)  # For data manipulation
library(ggplot2)     # For plotting (optional)
```

### Example 1: Response Spectra Calculation

```r
# Create sample acceleration data (2 Hz sine wave)
dt <- 0.01  # 100 Hz sampling rate
t <- seq(0, 10, by = dt)
acceleration <- 0.2 * sin(2 * pi * 2 * t) * exp(-0.1 * t)  # Damped oscillation

# Calculate response spectra
spectra <- get_Spectra(
    s = acceleration,
    t = t,
    xi = 0.05  # 5% damping
)

# Plot results
plot(spectra$Tn, spectra$PSA, 
     type = "l", log = "xy",
     xlab = "Period (s)", ylab = "PSA (m/s²)",
     main = "Pseudo-Spectral Acceleration")
grid(TRUE)
```

### Example 2: FFT Spectrum Analysis

```r
# Multi-frequency signal
signal <- sin(2 * pi * 5 * t) + 0.5 * sin(2 * pi * 15 * t) + 0.1 * rnorm(length(t))

# Compute FFT spectrum
spectrum <- build_FFT(
    s = signal,
    t = t,
    zp = 256,    # Zero padding
    kf = 0.8     # Frequency cutoff factor
)

# Plot frequency spectrum
plot(spectrum$fs, spectrum$A,
     type = "l",
     xlab = "Frequency (Hz)", ylab = "Amplitude",
     main = "FFT Amplitude Spectrum")
```

### Example 3: Empirical Mode Decomposition

```r
# Decompose signal into intrinsic mode functions
imfs <- get_imf(
    s = signal,
    ts = t,
    method = "eemd",     # Ensemble EMD
    trials = 20,         # Number of ensemble trials
    noise.amp = 0.1,     # Noise amplitude
    max.imf = 8          # Maximum number of IMFs
)

# Plot first few IMFs
library(ggplot2)
ggplot(imfs[IMF %in% c("IMF1", "IMF2", "IMF3", "signal")], 
       aes(x = t, y = s, color = IMF)) +
    geom_line() +
    facet_wrap(~IMF, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(title = "Empirical Mode Decomposition",
         x = "Time (s)", y = "Amplitude")
```

### Example 4: Intensity Measures Calculation

```r
# Prepare time series data in required format
ts_data <- data.table(
    RecordSN = 1,
    DIR = "H1",
    OCID = "EQ001",
    ID = "AT",  # Acceleration time series
    s = acceleration,
    t = t
)

# Calculate comprehensive intensity measures
intensity <- get_Intensity(ts_data, TargetUnits = "mm")

# View results
print(intensity$ITW)  # Wide format results
```

### Example 5: Newmark Displacement Analysis

```r
# Calculate sliding displacement for slope stability
kh <- 0.15  # Horizontal seismic coefficient
displacement <- get_ND(
    AT = acceleration * 9806.65,  # Convert to mm/s²
    t = t,
    kh = kh,
    TOL = 1e-3,
    full = TRUE  # Return full time history
)

# Plot displacement time history
plot(t, displacement * 1000,  # Convert to mm
     type = "l",
     xlab = "Time (s)", ylab = "Displacement (mm)",
     main = paste("Newmark Displacement (kh =", kh, ")"))
```

## Data Formats

The package is designed to work with both:

1. **Separate vectors**: `s` (signal) and `t` (time) parameters
2. **data.table format**: `.x` parameter with columns 's' and 't'

```r
# Method 1: Separate vectors
result1 <- get_Spectra(s = acceleration, t = time_vector)

# Method 2: data.table
acc_data <- data.table(s = acceleration, t = time_vector)
result2 <- get_Spectra(.x = acc_data)
```

## Advanced Features

### Spectral Smoothing
```r
# Smooth response spectra using different methods
smoothed <- smooth_Spectra(
    Tn = spectra$Tn,
    A = spectra$PSA,
    method = "sg",     # Savitzky-Golay filter
    window = 5,        # Window size
    po = 3            # Polynomial order
)
```

### IMF Filtering
```r
# Remove high-frequency IMFs (noise reduction)
filtered_signal <- remove_IMF(
    X = original_signal,
    EMD = emd_result,
    removeIMF1 = 2,    # Remove first 2 IMFs
    removeIMFn = 0     # Keep all remaining IMFs
)
```

## Function Reference

### Core Functions
- `get_Spectra()`: Calculate response spectra
- `build_FFT()`: Compute FFT spectrum
- `get_imf()`: Empirical mode decomposition
- `get_Intensity()`: Calculate intensity measures
- `get_ND()`: Newmark displacement analysis

### Utility Functions
- `smooth_Spectra()`: Spectral smoothing
- `remove_IMF()`: IMF filtering
- `map_OCID()`: OCID mapping utilities
- Various internal processing functions

## Applications

This package is particularly useful for:

- **Earthquake Engineering**: Response spectra calculation, intensity measures
- **Structural Dynamics**: Modal analysis, system identification
- **Geotechnical Engineering**: Slope stability analysis using Newmark method
- **Signal Processing**: Non-stationary signal analysis, time-frequency decomposition
- **Research**: Advanced signal processing for ground motion studies

## Performance Notes

- Functions are optimized for large datasets using `data.table`
- Matrix exponential calculations use efficient algorithms from the `expm` package
- Ensemble methods (EEMD, CEEMD) can be computationally intensive for long signals
- Memory cleanup is implemented using `on.exit()` for better memory management

## Contributing

We welcome contributions! Please see our contributing guidelines and feel free to submit issues or pull requests.

## Citation

If you use this package in your research, please cite:

```
Verri Kozlowski, A. (2024). gmsp: Ground Motion Signal Processing Tools. 
R package version 0.3.0.
```

## License

This project is licensed under the terms specified in the LICENSE file.

## Support

For questions, bug reports, or feature requests, please open an issue on the GitHub repository or contact the maintainer at averri@fi.uba.ar.

## Acknowledgments

This package builds upon several excellent R packages for signal processing and numerical computation, including `EMD`, `hht`, `signal`, `pracma`, and `data.table`. We thank the developers of these packages for their contributions to the R ecosystem.