# * ----
buildNavPanel.CEEMD <- function(id){
  ns <- NS(id)
  nav_panel(
    title="CEEMD", value=id,
    withMathJax(),
    navset_tab(
      nav_panel(title="ATH (H1)",navPanel.CEEMD.AT.H1(id)),
      nav_panel(title="ATH (H2)",navPanel.CEEMD.AT.H2(id)),
      nav_panel(title="ATH (UP)",navPanel.CEEMD.AT.UP(id)),
      nav_panel(title="VTH (H1)",navPanel.CEEMD.VT.H1(id)),
      nav_panel(title="VTH (H2)",navPanel.CEEMD.VT.H2(id)),
      nav_panel(title="VTH (UP)",navPanel.CEEMD.VT.UP(id)),
      nav_panel(title="DTH (H1)",navPanel.CEEMD.DT.H1(id)),
      nav_panel(title="DTH (H2)",navPanel.CEEMD.DT.H2(id)),
      nav_panel(title="DTH (UP)",navPanel.CEEMD.DT.UP(id)),
    )
  )# nav_panel
}

navPanel.CEEMD.AT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration (AT.H1)"),
    column(width = 12,highchartOutput(ns("CEEMD.AT.H1"),height = "700px"))
  )
}
navPanel.CEEMD.AT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration  (AT.H2)"),
    column(width = 12,highchartOutput(ns("CEEMD.AT.H2"),height = "700px"))
  )
}
navPanel.CEEMD.AT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Acceleration  (AT.UP)"),
    column(width = 12,highchartOutput(ns("CEEMD.AT.UP"),height = "700px"))
  )
}

navPanel.CEEMD.VT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H1)"),
    column(width = 12,highchartOutput(ns("CEEMD.VT.H1"),height = "700px"))
  )
}
navPanel.CEEMD.VT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H2)"),
    column(width = 12,highchartOutput(ns("CEEMD.VT.H2"),height = "700px"))
  )
}
navPanel.CEEMD.VT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Velocity (VT.UP)"),
    column(width = 12,highchartOutput(ns("CEEMD.VT.UP"),height = "700px"))
  )
}

navPanel.CEEMD.DT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H1)"),
    column(width = 12,highchartOutput(ns("CEEMD.DT.H1"),height = "700px"))
  )
}
navPanel.CEEMD.DT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H2)"),
    column(width = 12,highchartOutput(ns("CEEMD.DT.H2"),height = "700px"))
  )
}
navPanel.CEEMD.DT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Displacements (DT.UP)"),
    column(width = 12,highchartOutput(ns("CEEMD.DT.UP"),height = "700px"))
  )
}

