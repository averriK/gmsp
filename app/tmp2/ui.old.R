
ui <- fluidPage(
  # shinyjs::useShinyjs(),
  
  
  # withMathJax(),
  # titlePanel("GMRS"),
  
  page_navbar(
    header=NULL,
    footer=NULL,#"Ground Motion Record Selection",
    lang="es",
    theme = bs_theme(bootswatch = "united") , # bootswatch_themes(5)
    title = "GMRS",
    id = "nav",
    selected = "TS",
    sidebar = sidebar(
      position = "left",
      width = 300,
      open="desktop",
      fillable=TRUE,
      title = NULL,
      conditionalPanel("input.nav == 'TS'", id="TS", TS.sidebar(id="TS")),
      conditionalPanel("input.nav == 'FFT'", id="FFT", FFT.sidebar(id="FFT")),
      conditionalPanel("input.nav == 'PSA'", id="PSA", PSA.sidebar(id="PSA")),
      conditionalPanel("input.nav == 'STFT'", id="STFT", STFT.sidebar(id="STFT")),
      conditionalPanel("input.nav == 'WTC'", id="WTC", WTC.sidebar(id="WTC")),
      conditionalPanel("input.nav == 'EMD'", id="EMD", EMD.sidebar(id="EMD"))
    ),
    # Time Series Panel ----
    TS.ui(id="TS"),
    
    # Fourier Transform Panel ----
    FFT.ui(id="FFT"),
    
    # Response Spectra ----
    PSA.ui(id="PSA"),
    
    # Short-Time Fourier Transform Panel ----
    STFT.ui(id="STFT"),
    
    # Wavelet Transform Panel ----
    WTC.ui(id="WTC"),
    
    # EMD  Panel ----
    EMD.ui(id="EMD")
  )
)
