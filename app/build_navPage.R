
# -----------------------------------------------------------------------------
GMDB.navPage <- function(id) {
  tabPanel(
    "GMDB", 
    sidebarLayout(
      sidebarPanel(
        # GMDB specific sidebar elements
      ),
      mainPanel(
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
        TS.sidebarPanel("TS")
      ),
      mainPanel(
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
        SDOF.sidebarPanel("SDOF")
      ),
      mainPanel(
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
        STFT.sidebarPanel("STFT")
      ),
      mainPanel(
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
        FFT.sidebarPanel("FFT")
      ),
      mainPanel(
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
        CWT.sidebarPanel("CWT")
      ),
      mainPanel(
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
        EMD.sidebarPanel("EMD")
      ),
      mainPanel(
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
