# * ----
buildNavPanel.WTC <- function(id){
  ns <- NS(id)
  nav_panel(
    title="WTC", value=id,
    withMathJax(),
    navset_tab(
      nav_panel(title="ATH (H1)",navPanel.WTC.AT.H1(id)),
      nav_panel(title="ATH (H2)",navPanel.WTC.AT.H2(id)),
      nav_panel(title="ATH (UP)",navPanel.WTC.AT.UP(id)),
      nav_panel(title="VTH (H1)",navPanel.WTC.VT.H1(id)),
      nav_panel(title="VTH (H2)",navPanel.WTC.VT.H2(id)),
      nav_panel(title="VTH (UP)",navPanel.WTC.VT.UP(id)),
      nav_panel(title="DTH (H1)",navPanel.WTC.DT.H1(id)),
      nav_panel(title="DTH (H2)",navPanel.WTC.DT.H2(id)),
      nav_panel(title="DTH (UP)",navPanel.WTC.DT.UP(id)),
    )
  )# nav_panel
}

navPanel.WTC.AT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration (AT.H1)"),
    column(width = 12,highchartOutput(ns("WTC.AT.H1"),height = "700px"))
  )
}
navPanel.WTC.AT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration  (AT.H2)"),
    column(width = 12,highchartOutput(ns("WTC.AT.H2"),height = "700px"))
  )
}
navPanel.WTC.AT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Acceleration  (AT.UP)"),
    column(width = 12,highchartOutput(ns("WTC.AT.UP"),height = "700px"))
  )
}

navPanel.WTC.VT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H1)"),
    column(width = 12,highchartOutput(ns("WTC.VT.H1"),height = "700px"))
  )
}
navPanel.WTC.VT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H2)"),
    column(width = 12,highchartOutput(ns("WTC.VT.H2"),height = "700px"))
  )
}
navPanel.WTC.VT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Velocity (VT.UP)"),
    column(width = 12,highchartOutput(ns("WTC.VT.UP"),height = "700px"))
  )
}

navPanel.WTC.DT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H1)"),
    column(width = 12,highchartOutput(ns("WTC.DT.H1"),height = "700px"))
  )
}
navPanel.WTC.DT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H2)"),
    column(width = 12,highchartOutput(ns("WTC.DT.H2"),height = "700px"))
  )
}
navPanel.WTC.DT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Displacements (DT.UP)"),
    column(width = 12,highchartOutput(ns("WTC.DT.UP"),height = "700px"))
  )
}

