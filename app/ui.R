source("setup.R",local = TRUE)
source("buildUI.R")
ui <- navbarPage(
  "Shiny Dashboard",
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
  ),
  # GMDB.tabPanel("GMDB"),
  TS.tabPanel("TS"),
  SDOF.tabPanel("SDOF"),
  FFT.tabPanel("FFT"),
  STFT.tabPanel("STFT"),
  CWT.tabPanel("CWT"),
  EMD.tabPanel("EMD"),
  
)


