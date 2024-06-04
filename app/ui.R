library(shiny)
library(shinydashboard)


# Define the UI
ui <- dashboardPage(
  header = dashboardHeader(),
  sidebar = NULL,
  
  body = dashboardBody(
    TS.tabPanel(),
    SDOF.tabPanel()
    
  ), # dashboardBody
  
  title = "Grount-Motion Signal Processing (GMSP)"
)

