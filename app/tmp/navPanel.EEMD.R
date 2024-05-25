# * ----
buildNavPanel.EEMD <- function(id){
  ns <- NS(id)
  nav_panel(
    title="EEMD", value=id,
    withMathJax(),
    navset_tab(
      nav_panel(title="ATH (H1)",navPanel.EEMD.AT.H1(id)),
      nav_panel(title="ATH (H2)",navPanel.EEMD.AT.H2(id)),
      nav_panel(title="ATH (UP)",navPanel.EEMD.AT.UP(id)),
      nav_panel(title="VTH (H1)",navPanel.EEMD.VT.H1(id)),
      nav_panel(title="VTH (H2)",navPanel.EEMD.VT.H2(id)),
      nav_panel(title="VTH (UP)",navPanel.EEMD.VT.UP(id)),
      nav_panel(title="DTH (H1)",navPanel.EEMD.DT.H1(id)),
      nav_panel(title="DTH (H2)",navPanel.EEMD.DT.H2(id)),
      nav_panel(title="DTH (UP)",navPanel.EEMD.DT.UP(id)),
    )
  )# nav_panel
}

navPanel.EEMD.AT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration (AT.H1)"),
    column(width = 12,highchartOutput(ns("EEMD.AT.H1"),height = "700px"))
  )
}
navPanel.EEMD.AT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration  (AT.H2)"),
    column(width = 12,highchartOutput(ns("EEMD.AT.H2"),height = "700px"))
  )
}
navPanel.EEMD.AT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Acceleration  (AT.UP)"),
    column(width = 12,highchartOutput(ns("EEMD.AT.UP"),height = "700px"))
  )
}

navPanel.EEMD.VT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H1)"),
    column(width = 12,highchartOutput(ns("EEMD.VT.H1"),height = "700px"))
  )
}
navPanel.EEMD.VT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H2)"),
    column(width = 12,highchartOutput(ns("EEMD.VT.H2"),height = "700px"))
  )
}
navPanel.EEMD.VT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Velocity (VT.UP)"),
    column(width = 12,highchartOutput(ns("EEMD.VT.UP"),height = "700px"))
  )
}

navPanel.EEMD.DT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H1)"),
    column(width = 12,highchartOutput(ns("EEMD.DT.H1"),height = "700px"))
  )
}
navPanel.EEMD.DT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H2)"),
    column(width = 12,highchartOutput(ns("EEMD.DT.H2"),height = "700px"))
  )
}
navPanel.EEMD.DT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Displacements (DT.UP)"),
    column(width = 12,highchartOutput(ns("EEMD.DT.UP"),height = "700px"))
  )
}

