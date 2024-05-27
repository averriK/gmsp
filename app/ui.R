library(shiny)
library(shinydashboard)


# Define UI module for Navigation Panel
navbarPageUI <- function(id,title="Time Series (TS)",hc=TRUE,title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  if(hc==TRUE) {
    UI <- navbarPage(
      title=title,
      id = id,
      fluidRow(
        column(
          width=8,
          tabsetPanel(
            tabPanel(title=title.AT,highchartOutput(ns("AT"),height = height)),
            tabPanel(title=title.VT,highchartOutput(ns("VT"),height = height)),
            tabPanel(title=title.DT,highchartOutput(ns("DT"),height = height))
          ) # tabsetPanel
        ) # column
      ) # fluidRow
    ) # navbarPage
  }  
  
  if(hc==FALSE){
    UI <-navbarPage(
      title=title,
      id = id,
      fluidRow(
        column(
          width=8,
          tabsetPanel(
            tabPanel(title=title.AT,plotOutput(ns("AT"),height = height)),
            tabPanel(title=title.VT,plotOutput(ns("VT"),height = height)),
            tabPanel(title=title.DT,plotOutput(ns("DT"),height = height))
          ) # tabsetPanel
        ) # column
      ) # fluidRow
    ) # navbarPage
  }
  return(UI)
  
}


# Define the UI
ui <- dashboardPage(
  header = dashboardHeader(),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("GMDB", tabName = "gmdb"),
      menuItem("TS", tabName = "ts"),
      menuItem("FFT", tabName = "fft"),
      menuItem("STFT", tabName = "stft"),
      menuItem("WTC", tabName = "wtc"),
      menuItem("SDOF", tabName = "sdof"),
      menuItem("EMD", tabName = "emd")
    )
  ),
  
  body = dashboardBody(
    tabItems(
      # tabItem(tabName = "gmdb"), # tabItem
      tabItem(tabName = "ts",
              navbarPageUI(
                id="TS",hc=TRUE,title="Time Series (TS)",
                title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)")), # tabItem
      tabItem(tabName = "fft",
              navbarPageUI(
                id="FFT",hc=TRUE,title="Fast Fourier Transform (FFT)",
                title.AT="FFT(AT)",title.VT="FFT(VT)",title.DT="FFT(DT)")), # tabItem
      tabItem(tabName = "stft",
              navbarPageUI(
                id="STFT",hc=FALSE,title="Short-Time Fourier Transform (STFT)",
                title.AT="STFT(AT)",title.VT="STFT(VT)",title.DT="STFT(DT)")), # tabItem
      tabItem(tabName = "wtc",
              navbarPageUI(
                id="WTC",hc=FALSE,title="Continuous Wavelet Transform (CWT)",
                title.AT="CWT(AT)",title.VT="CWT(VT)",title.DT="CWT(DT)")), # tabItem
      tabItem(tabName = "sdof",
              navbarPageUI(
                id="SDOF",hc=TRUE,title="Single Degree of Freedom (SDOF)",
                title.AT="PSA",title.VT="PSV",title.DT="SD")), # tabItem
      tabItem(tabName = "emd",
              navbarPageUI(
                id="EMD",hc=TRUE,title="Empirical Mode Decomposition (EMD)",
                title.AT="EMD(AT)",title.VT="EMD(VT)",title.DT="EMD(DT)")) #tabItem
      
    ) # tabItems
  ), # dashboardBody
  
  title = "ShinyDashboard Example"
)

