# * ----
buildNavPanel.EMD <- function(id){
  ns <- NS(id)
  nav_panel(
    title="EMD", value=id,
    withMathJax(),
    navset_tab(
      nav_panel(title="ATH (H1)",navPanel.EMD.AT.H1(id)),
      nav_panel(title="ATH (H2)",navPanel.EMD.AT.H2(id)),
      nav_panel(title="ATH (UP)",navPanel.EMD.AT.UP(id)),
      nav_panel(title="VTH (H1)",navPanel.EMD.VT.H1(id)),
      nav_panel(title="VTH (H2)",navPanel.EMD.VT.H2(id)),
      nav_panel(title="VTH (UP)",navPanel.EMD.VT.UP(id)),
      nav_panel(title="DTH (H1)",navPanel.EMD.DT.H1(id)),
      nav_panel(title="DTH (H2)",navPanel.EMD.DT.H2(id)),
      nav_panel(title="DTH (UP)",navPanel.EMD.DT.UP(id)),
    )
  )# nav_panel
}

navPanel.EMD.AT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration (AT.H1)"),
    column(width = 12,highchartOutput(ns("EMD.AT.H1"),height = "700px"))
  )
}
navPanel.EMD.AT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration  (AT.H2)"),
    column(width = 12,highchartOutput(ns("EMD.AT.H2"),height = "700px"))
  )
}
navPanel.EMD.AT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Acceleration  (AT.UP)"),
    column(width = 12,highchartOutput(ns("EMD.AT.UP"),height = "700px"))
  )
}

navPanel.EMD.VT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H1)"),
    column(width = 12,highchartOutput(ns("EMD.VT.H1"),height = "700px"))
  )
}
navPanel.EMD.VT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H2)"),
    column(width = 12,highchartOutput(ns("EMD.VT.H2"),height = "700px"))
  )
}
navPanel.EMD.VT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Velocity (VT.UP)"),
    column(width = 12,highchartOutput(ns("EMD.VT.UP"),height = "700px"))
  )
}

navPanel.EMD.DT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H1)"),
    column(width = 12,highchartOutput(ns("EMD.DT.H1"),height = "700px"))
  )
}
navPanel.EMD.DT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H2)"),
    column(width = 12,highchartOutput(ns("EMD.DT.H2"),height = "700px"))
  )
}
navPanel.EMD.DT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Displacements (DT.UP)"),
    column(width = 12,highchartOutput(ns("EMD.DT.UP"),height = "700px"))
  )
}

