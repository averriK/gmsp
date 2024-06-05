source("setup.R",local = TRUE)
source("build_sidebarPanel.R")
source("build_mainPanel.R")
source("build_navPage.R")
ui <- navbarPage(
  title="GMRS",
  header = NULL,
  
  GMDB.navPage("GMDB"),
  TS.navPage("TS"),
  SDOF.navPage("SDOF"),
  FFT.navPage("FFT"),
  STFT.navPage("STFT"),
  CWT.navPage("CWT"),
  EMD.navPage("EMD"),
  
)


