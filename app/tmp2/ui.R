library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(htmltools)

# Define the UI
ui <- navbarPage(
  
  TS.tabPanel("TS"),
  SDOF.tabPanel("SDOF"),
  FFT.tabPanel("FFT"),
  STFT.tabPanel("STFT"),
  CWT.tabPanel("CWT"),
  EMD.tabPanel("EMD")
  
)

