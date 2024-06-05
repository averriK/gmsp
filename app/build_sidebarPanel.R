TS.sidebarPanel <- function(id) {
  ns <- NS(id)
  tagList(
    helpText("Direction:"),
    prettyRadioButtons(
      inputId = ns("ocid"),
      label = NULL,
      choices = c("H1", "H2","UP"),
      selected = "H1",
      inline = TRUE,
      status = "danger",
      fill = TRUE
    )
  )
}
STFT.sidebarPanel <- function(id){
  ns=NS(id)
  tagList(
    helpText("Windows Length k=N/NW:"),
    sliderTextInput(
      inputId = ns("k"),
      label = NULL, 
      choices = c(4,8,16,32,64,128,256,512),
      selected = 32,
      grid = TRUE),
    
    
    # Zero Padding----
    helpText("Zero padding zp:"),
    sliderTextInput(
      inputId = ns("zp"),
      label = NULL, 
      choices = c(128,256,512,1024,2048,4096,8192),
      selected = 256,
      grid = TRUE),
    
    # Windows Overlap ----
    
    helpText("Windows Overlap [%]:"),
    sliderTextInput(
      inputId = ns("ovlp"),
      label = NULL, 
      choices = c(10,25,50,75,90),
      selected = 75,
      grid = TRUE),
    
    # Nyqist ----
    
    helpText("Max. Freq. kf=fmax/fNyq:)"),
    sliderTextInput(
      inputId = ns("kf"),
      label = NULL, 
      choices = c(seq(0.01,0.14,by=0.01),seq(0.15,0.95,by=0.05)),
      selected = 0.25,
      grid = FALSE),
    
    # norm ----
    helpText("Normalize:"),
    radioButtons(
      inputId = ns("norm"),
      label = NULL,
      choices = c(TRUE,FALSE),
      selected = TRUE,
      inline = TRUE
    ),
    
    # osc ----
    helpText("Oscillator:"),
    radioButtons(
      inputId = ns("osc"),
      label = NULL,
      choices = c(TRUE,FALSE),
      selected = TRUE,
      inline = TRUE
    ),
    
    # scale ----
    helpText("Scale:"),
    radioButtons(
      inputId = ns("scale"),
      label = NULL,
      choices = c(TRUE,FALSE),
      selected = TRUE,
      inline = TRUE
    ),
    # OCID ----
    helpText("OCID:"),
    radioButtons(
      inputId = ns("ocid"),
      label = NULL,
      choices = c("H1", "H2","UP"),
      selected = "H1",
      inline = TRUE
    )
    
  ) # tagList
  
  
} # STFT.sidebar
SDOF.sidebarPanel <- function(id = "SDOF"){
  ns=NS(id)
  tagList(
    helpText("Smoothing Method:"),
    prettyRadioButtons(
      inputId = ns("smoothing"),
      label = NULL,
      choiceValues  = c("none", "ma","sg","ema","sm","gk"),
      choiceNames = c("No Smoothing", "Moving Average","Savitzky-Golay","Exponential Moving Average","Spline Smooth","Gasussian Kernel"),
      selected = "none",
      inline = FALSE,
      status = "danger",
      fill = TRUE
    ),
    
    helpText("Window Size (smooth level):"),
    sliderTextInput(
      inputId = ns("window"),
      label = NULL, 
      choices = seq(2,32),
      selected = 8,
      grid = TRUE),
    
    helpText("Cross-Validation (spline):"),
    prettyRadioButtons(
      inputId = ns("cv"),
      label = NULL,
      choices = c(TRUE,FALSE),
      selected = FALSE,
      inline = TRUE,
      status = "danger",
      fill = TRUE
    ),
    
    
    helpText("Sigma (Gauss Kernel):"),
    sliderTextInput(
      inputId = ns("sigma"),
      label = NULL, 
      choices = seq(0.5,5,by=0.1),
      selected = 2,
      grid = TRUE),
    
    helpText("Critical Damping Ratio:"),
    sliderTextInput(
      inputId = ns("xi"),
      label = NULL, 
      choices = seq(0.01,0.30,by=0.01),
      selected = 0.05,
      grid = TRUE),
    
    helpText("OCID:"),
    prettyRadioButtons(
      inputId = ns("ocid"),
      label = NULL,
      choices = c("H1", "H2","UP"),
      selected = "H1",
      inline = TRUE,
      status = "danger",
      fill = TRUE
    )
  )
}
FFT.sidebarPanel <- function(id){
  ns=NS(id)
  tagList(
    
    helpText("Smoothing Method:"),
    prettyRadioButtons(
      inputId = ns("smoothing"),
      label = NULL,
      choiceValues  = c("none", "ma","sg","ema","sm","gk"),
      choiceNames = c("No Smoothing", "Moving Average","Savitzky-Golay","Exponential Moving Average","Spline Smooth","Gasussian Kernel"),
      selected = "none",
      inline = FALSE,
      status = "danger",
      fill = TRUE
    ),
    
    helpText("Window Size (smooth level):"),
    sliderTextInput(
      inputId = ns("window"),
      label = NULL, 
      choices = seq(2,32),
      selected = 8,
      grid = TRUE),
    
    helpText("Cross-Validation (spline):"),
    prettyRadioButtons(
      inputId = ns("cv"),
      label = NULL,
      choices = c(TRUE,FALSE),
      selected = FALSE,
      inline = TRUE,
      status = "danger",
      fill = TRUE
    ),
    
    
    helpText("Sigma (Gauss Kernel):"),
    sliderTextInput(
      inputId = ns("sigma"),
      label = NULL, 
      choices = seq(0.5,5,by=0.1),
      selected = 2,
      grid = TRUE),
    
    # Nyqist ----
    
    helpText("Max. Freq. kf=fmax/fNyq:)"),
    sliderTextInput(
      inputId = ns("kf"),
      label = NULL, 
      choices = c(seq(0.01,0.14,by=0.01),seq(0.15,0.95,by=0.05)),
      selected = 0.25,
      grid = FALSE),
    
    # Zero Padding----
    helpText("Zero padding zp:"),
    sliderTextInput(
      inputId = ns("zp"),
      label = NULL, 
      choices = c(64,128,256,512,1024,2048,4096),
      selected = 256,
      grid = TRUE),
    
    
    
    helpText("OCID:"),
    prettyRadioButtons(
      inputId = ns("ocid"),
      label = NULL,
      choices = c("H1", "H2","UP"),
      selected = "H1",
      inline = TRUE,
      status = "danger",
      fill = TRUE
    )
  )
}
EMD.sidebarPanel <- function(id){
  ns=NS(id)
  tagList(
    helpText("Max number of IMFs:"),
    sliderTextInput(
      inputId = ns("max.imf"),
      label = NULL, 
      choices = seq(1,20),
      selected = 12,
      grid = TRUE
    ),
    
    
    
    helpText("Boundary type:"),
    prettyRadioButtons(
      inputId = ns("boundary"),
      label = NULL,
      choices=c("none","wave","symmetric","periodic","evenodd"),
      selected="wave",
      inline = FALSE,
      status = "danger",
      fill = TRUE
    ),
    
    helpText("Stop rule:"),
    prettyRadioButtons(
      inputId = ns("stop.rule"),
      label = NULL,
      choices=c("type1","type2","type3","type4","type5"),
      selected="type5",
      inline = FALSE,
      status = "danger",
      fill = TRUE
    ),
    
    
    # OCID ----
    helpText("OCID:"),
    prettyRadioButtons(
      inputId = ns("ocid"),
      label = NULL,
      choices = c("H1", "H2","UP"),
      selected = "H1",
      inline = TRUE,
      status = "danger",
      fill = TRUE
    )
    
  ) # tagList
  
  
} # EMD.sidebar


CWT.sidebarPanel <- function(id){
  ns=NS(id)
  tagList(
    helpText("Wavelet family:"),
    prettyRadioButtons(
      inputId = ns("mother"),
      label = NULL,
      choices=c("morlet","dog","paul"),
      selected="morlet",
      inline = TRUE, 
      status = "danger",
      fill = TRUE
    ),
    
    
    # Wavelet Scale ----
    helpText("Wavelet Scale j=1/dj:"),
    sliderTextInput(
      inputId = ns("j"),
      label = NULL, 
      choices = c(2,4,8,12,16,24,32,48,64,72,128,256,512),
      selected = 64,
      grid = TRUE),
    
    
    
    # norm ----
    helpText("Include Phases:"),
    prettyRadioButtons(
      inputId = ns("phase"),
      label = NULL,
      choices = c(TRUE,FALSE),
      selected = TRUE,
      inline = TRUE,
      status = "danger",
      fill = TRUE
    ),
    
    
    # OCID ----
    helpText("OCID:"),
    prettyRadioButtons(
      inputId = ns("ocid"),
      label = NULL,
      choices = c("H1", "H2","UP"),
      selected = "H1",
      inline = TRUE,
      status = "danger",
      fill = TRUE
    )
    
  ) # tagList
  
  
} # CWT.sidebar
