# * ----
buildNavPanel.STFT <- function(id){
  ns <- NS(id)
  nav_panel(
    title="STFT", value=id,
    withMathJax(),
    navset_tab(
      nav_panel(title="ATH (H1)",navPanel.STFT.AT.H1(id)),
      nav_panel(title="ATH (H2)",navPanel.STFT.AT.H2(id)),
      nav_panel(title="ATH (UP)",navPanel.STFT.AT.UP(id)),
      nav_panel(title="VTH (H1)",navPanel.STFT.VT.H1(id)),
      nav_panel(title="VTH (H2)",navPanel.STFT.VT.H2(id)),
      nav_panel(title="VTH (UP)",navPanel.STFT.VT.UP(id)),
      nav_panel(title="DTH (H1)",navPanel.STFT.DT.H1(id)),
      nav_panel(title="DTH (H2)",navPanel.STFT.DT.H2(id)),
      nav_panel(title="DTH (UP)",navPanel.STFT.DT.UP(id)),
    )
  )# nav_panel
}

navPanel.STFT.AT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration (AT.H1)"),
    column(width = 12,highchartOutput(ns("STFT.AT.H1"),height = "700px"))
  )
}
navPanel.STFT.AT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Acceleration  (AT.H2)"),
    column(width = 12,highchartOutput(ns("STFT.AT.H2"),height = "700px"))
  )
}
navPanel.STFT.AT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Acceleration  (AT.UP)"),
    column(width = 12,highchartOutput(ns("STFT.AT.UP"),height = "700px"))
  )
}

navPanel.STFT.VT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H1)"),
    column(width = 12,highchartOutput(ns("STFT.VT.H1"),height = "700px"))
  )
}
navPanel.STFT.VT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Velocity (VT.H2)"),
    column(width = 12,highchartOutput(ns("STFT.VT.H2"),height = "700px"))
  )
}
navPanel.STFT.VT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Velocity (VT.UP)"),
    column(width = 12,highchartOutput(ns("STFT.VT.UP"),height = "700px"))
  )
}

navPanel.STFT.DT.H1 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H1)"),
    column(width = 12,highchartOutput(ns("STFT.DT.H1"),height = "700px"))
  )
}
navPanel.STFT.DT.H2 <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Horizontal Displacements (DT.H2)"),
    column(width = 12,highchartOutput(ns("STFT.DT.H2"),height = "700px"))
  )
}
navPanel.STFT.DT.UP <- function(id){
  ns <- NS(id)
  fluidRow(
    helpText("Vertical Displacements (DT.UP)"),
    column(width = 12,highchartOutput(ns("STFT.DT.UP"),height = "700px"))
  )
}

