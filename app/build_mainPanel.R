# -----------------------------------------------------------------------------
TS.navbarPageUI <- function(id,title="Time Series (TS)",title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  
    navbarPage(
      title=title,
      id = id,
      fluidRow(
        column(
          width=12,
          tabsetPanel(
            tabPanel(title=title.AT,highchartOutput(ns("AT"),height = height)),
            tabPanel(title=title.VT,highchartOutput(ns("VT"),height = height)),
            tabPanel(title=title.DT,highchartOutput(ns("DT"),height = height))
          ) # tabsetPanel
        ) # column
      ) # fluidRow
    ) # navbarPage

}

SDOF.navbarPageUI <- function(id,title="Single Degree of Freedom (SDOF)",title.AT="Pseudo Spectral Acceleration (PSA)",title.VT="Pseudo Spectral Velocity (PSV)",title.DT="Spectrl Displacement (SD)",height="500px"){
  ns <- NS(id)
  navbarPage(
    title=title,
    id = id,
    fluidRow(
      column(
        width=12,
        tabsetPanel(
          tabPanel(title=title.AT,highchartOutput(ns("AT"),height = height)),
          tabPanel(title=title.VT,highchartOutput(ns("VT"),height = height)),
          tabPanel(title=title.DT,highchartOutput(ns("DT"),height = height))
        )
      )
    )
  ) # navbarPage
}

STFT.navbarPageUI <- function(id,title="Short-Time Fourier Transform (STFT)",title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  navbarPage(
    title=title,
    id = id,
    fluidRow(
      column(
        width=12,
        tabsetPanel(
          tabPanel(title=title.AT,plotOutput(ns("AT"),height = height)),
          tabPanel(title=title.VT,plotOutput(ns("VT"),height = height)),
          tabPanel(title=title.DT,plotOutput(ns("DT"),height = height))
        )
      )
    )
  ) # navbarPage
}

FFT.navbarPageUI <- function(id,title="Fast Fourier Transform (FFT)",title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  navbarPage(
    title=title,
    id = id,
    fluidRow(
      column(
        width=12,
        tabsetPanel(
          tabPanel(title=title.AT,plotOutput(ns("AT"),height = height)),
          tabPanel(title=title.VT,plotOutput(ns("VT"),height = height)),
          tabPanel(title=title.DT,plotOutput(ns("DT"),height = height))
        )
      )
    )
  ) # navbarPage
}

CWT.navbarPageUI <- function(id,title="Continuous Wavelet Transform (CWT)",title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  navbarPage(
    title=title,
    id = id,
    fluidRow(
      column(
        width=12,
        tabsetPanel(
          tabPanel(title=title.AT,plotOutput(ns("AT"),height = height)),
          tabPanel(title=title.VT,plotOutput(ns("VT"),height = height)),
          tabPanel(title=title.DT,plotOutput(ns("DT"),height = height))
        )
      )
    )
  ) # navbarPage
}

EMD.navbarPageUI <- function(id,title="Empirical Mode Decomposition (EMD)",title.AT="Acceleration (AT)",title.VT="Velocity (VT)",title.DT="Displacement (DT)",height="500px"){
  ns <- NS(id)
  navbarPage(
    title=title,
    id = id,
    fluidRow(
      column(
        width=12,
        tabsetPanel(
          tabPanel(title=title.AT,plotOutput(ns("AT"),height = height)),
          tabPanel(title=title.VT,plotOutput(ns("VT"),height = height)),
          tabPanel(title=title.DT,plotOutput(ns("DT"),height = height))
        )
      )
    )
  ) # navbarPage
}
