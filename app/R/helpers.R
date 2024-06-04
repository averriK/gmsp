# Define UI module for Navigation Panel
buildTabPanel <- function(id) {
  tabPanel(
    title=id,
    sidebarLayout(
      sidebarPanel(
        TS.sidebarPanel(id=id)
      ),
      mainPanel(
        TS.mainPanel(id = id)
      )
    )
  )
}


buildMainPanel <- function(id,title="Time Series (TS)",render="highchartOutput",title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  if(render=="highchartOutput") {
    UI <- navbarPage(title=title,
                     id = id,
                     fluidRow(
                       column(
                         width=8,
                         tabsetPanel(
                           tabPanel(title=title.AT,highchartOutput(ns("AT"),height = height)),
                           tabPanel(title=title.VT,highchartOutput(ns("VT"),height = height)),
                           tabPanel(title=title.DT,highchartOutput(ns("DT"),height = height))
                         )
                       )
                     )
    ) # navbarPage
  }  
  
  if(render=="plotOutput"){
    UI <-navbarPage(title=title,
                    id = id,
                    fluidRow(
                      column(
                        width=8,
                        tabsetPanel(
                          tabPanel(title=title.AT,plotOutput(ns("AT"),height = height)),
                          tabPanel(title=title.VT,plotOutput(ns("VT"),height = height)),
                          tabPanel(title=title.DT,plotOutput(ns("DT"),height = height))
                        )
                      )
                    )
    ) # navbarPage
  }
  return(UI)
  
}
