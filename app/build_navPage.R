
# -----------------------------------------------------------------------------
GMDB.navPage <- function(id) {
  nav_panel(
    title = "GMDB", 
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          # GMDB specific sidebar elements
        ),
        mainPanel(
          width = 6,
          # GMDB specific main content
        )
      )
    )
  )
}
TS.navPage <- function(id) {
  nav_panel(
    title = "TS",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          TS.sidebarPanel("TS")
        ),
        mainPanel(
          width = 7,
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
  )
}
SDOF.navPage <- function(id) {
  nav_panel(
    title ="SDOF",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          SDOF.sidebarPanel("SDOF")
        ),
        mainPanel(
          width = 7,
          SDOF.navbarPageUI(
            id = "SDOF", 
            title = "Single Degree of Freedom (SDOF)",
            title.AT = "Acceleration (AT)", 
            title.VT = "Velocity (VT)", 
            title.DT = "Displacement (DT)"
          )
        )
      )
    )
  )
}
STFT.navPage <- function(id) {
  nav_panel(
    title="STFT",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          STFT.sidebarPanel("STFT")
        ),
        mainPanel(
          width = 7,
          STFT.navbarPageUI(
            id = "STFT", 
            title = "Short Time Fourier Transform (STFT)",
            title.AT = "Acceleration (AT)", 
            title.VT = "Velocity (VT)", 
            title.DT = "Displacement (DT)"
          )
        )
      )
    )
  )
}
FFT.navPage <- function(id) {
  nav_panel(
    title="FFT",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          FFT.sidebarPanel("FFT")
        ),
        mainPanel(
          width = 7,
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
  )
}
CWT.navPage <- function(id) {
  nav_panel(
    title="CWT",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          CWT.sidebarPanel("CWT")
        ),
        mainPanel(
          width = 7,
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
    
  )
}
EMD.navPage <- function(id) {
  nav_panel(
    title="EMD",
    fluidPage(
      sidebarLayout(
        sidebarPanel(
          width = 3,
          EMD.sidebarPanel("EMD")
        ),
        mainPanel(
          width = 7,
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
  )
}
