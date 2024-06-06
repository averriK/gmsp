
# -----------------------------------------------------------------------------
GMDB.navPage <- function(id) {
  tabPanel(
    "GMDB", 
    sidebarLayout(
      sidebarPanel(
        width=2,
        # GMDB specific sidebar elements
      ),
      mainPanel(
        width=7,
        # GMDB specific main content
      )
    )
  )
}

TS.navPage <- function(id) {
  tabPanel(
    "TS",
    sidebarLayout(
      sidebarPanel(
        width=2,
        TS.sidebarPanel("TS")
      ),
      mainPanel(
        width=7,
        TS.navbarPageUI(
          id = "TS", 
          title = "Time Series (TS)",
          title.AT = "Acceleration (AT)", 
          title.VT = "Velocity (VT)", 
          title.DT = "Displacement (DT)"
        )
      )
    )
  )
}
SDOF.navPage <- function(id) {
  tabPanel(
    "SDOF",
    sidebarLayout(
      sidebarPanel(
        width=2,
        SDOF.sidebarPanel("SDOF")
      ),
      mainPanel(
        width=7,
        SDOF.navbarPageUI(
          id = "SDOF", 
          title = "Single Degree of Freedom (SDOF)",
          title.AT = "Pseudo Spectral Acceleration (PSA)", 
          title.VT = "Pseudo Spectral Velocity (PSV)", 
          title.DT = "Spectral Displacement (SD)"
        )
      )
    )
  )
}
STFT.navPage <- function(id) {
  tabPanel(
    "STFT",
    sidebarLayout(
      sidebarPanel(
        width=2,
        STFT.sidebarPanel("STFT")
      ),
      mainPanel(
        width=7,
        STFT.navbarPageUI(
          id = "STFT", 
          title = "Short-Time Fourier Transform (STFT)",
          title.AT = "Acceleration (AT)", 
          title.VT = "Velocity (VT)", 
          title.DT = "Displacement (DT)"
        )
      )
    )
  )
}
FFT.navPage <- function(id) {
  tabPanel(
    "FFT",
    sidebarLayout(
      sidebarPanel(
        width=2,
        FFT.sidebarPanel("FFT")
      ),
      mainPanel(
        width=7,
        FFT.navbarPageUI(
          id = "FFT", 
          title = "Fast Fourier Transform (FFT)",
          title.AT = "Acceleration (AT)", 
          title.VT = "Velocity (VT)", 
          title.DT = "Displacement (DT)"
        )
      )
    )
  )
}
CWT.navPage <- function(id) {
  tabPanel(
    "CWT",
    sidebarLayout(
      sidebarPanel(
        width=2,
        CWT.sidebarPanel("CWT")
      ),
      mainPanel(
        width=7,
        CWT.navbarPageUI(
          id = "CWT", 
          title = "Continuous Wavelet Transform (CWT)",
          title.AT = "Acceleration (AT)", 
          title.VT = "Velocity (VT)", 
          title.DT = "Displacement (DT)"
        )
      )
    )
  )
}
EMD.navPage <- function(id) {
  tabPanel(
    "EMD",
    sidebarLayout(
      sidebarPanel(
        width=2,
        EMD.sidebarPanel("EMD")
      ),
      mainPanel(
        width=7,
        EMD.navbarPageUI(
          id = "EMD", 
          title = "Empirical Mode Decomposition (EMD)",
          title.AT = "Acceleration (AT)", 
          title.VT = "Velocity (VT)", 
          title.DT = "Displacement (DT)"
        )
      )
    )
  )
}
