source("setup.R",local = TRUE)
source("build_sidebarPanel.R")
source("build_mainPanel.R")
source("build_navPage.R")
ui <- navbarPage(
  title="GMRS",
  header = tagList(
    # GMRS.header()
  ),
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


