source("setup.R",local = TRUE)
source("global.R",local = TRUE)
source("helpers.R",local = TRUE)

source("R/GMDB.R",local = TRUE)
source("R/TS.R",local = TRUE)
source("R/SDOF.R",local = TRUE)
source("R/FFT.R",local = TRUE)
source("R/STFT.R",local = TRUE)
source("R/WTC.R",local = TRUE)
source("R/EMD.R",local = TRUE)


server <- function(input, output, session) {
  output$mainContent <- renderUI({
    switch(
      input$topmenu,
      "ts" = navbarPageUI(
        id = "TS", hc = TRUE, title = "Time Series (TS)",
        title.AT = "Acceleration (AT)", title.VT = "Velocity (VT)", title.DT = "Displacement (DT)"
      ),
      "sdof" = SDOF.navbarPageUI(
        id = "SDOF", title = "Single Degree of Freedom (SDOF)",
        title.AT = "Pseudo Spectral Acceleration (PSA)", title.VT = "Pseudo Spectral Velocity (PSV)", title.DT = "Spectral Displacement (SD)"
      )
      # Add further cases here if needed
    )
  })
  
  TS.server(id = "TS", .data = TSL, series = "AT", yAxis.legend = "A(t)")
  TS.server(id = "TS", .data = TSL, series = "VT", yAxis.legend = "V(t)")
  TS.server(id = "TS", .data = TSL, series = "DT", yAxis.legend = "D(t)")
  
  SDOF.server(id = "SDOF", .data = TSL)
}


