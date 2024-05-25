STFT.sidebar <- function(id){
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

STFT.server <- function(id,.data,series){
  moduleServer(
    id, 
    function(input, output, session){
      stopifnot(series %in% c("AT","VT","DT"))
      observe({
        cat("Current OCID value:", input$ocid, "\n")  
        cat("Current kf value:", input$kf, "\n")  
        cat("Current k value:", input$k, "\n")  
        cat("Current ovlp value:", input$ovlp, "\n")  
      })
      
      DATA <- reactive({
        req(input$ocid)
        stopifnot(length(input$ocid)==1 && input$ocid %in% c("H1","H2","UP"))
        req(.data)
        DT <- .data[ID %in% series & OCID %in% input$ocid,.(t,s)]
        if(nrow(DT)==0) return(NULL)
        k <- req(input$k) |> as.integer()
        wl <- nrow(DT)/k
        while (wl %% 2 != 0) {
          DT <- DT[-1]  # Remove the first row from DT
          wl <- nrow(DT)/k
        }
        
        DT
        
      })
      
      WAVE <- reactive({
        req(DATA())
        k <- req(input$k) |> as.integer()
        wl <- nrow(DATA())/k
        dt <- mean(diff(DATA()$t))
        DT <- tuneR::Wave(left=DATA()$s, samp.rate=1/dt,bit=32)
        DT
      })
      
    
  
      
      output[[series]] <- renderPlot({
        req(DATA())
        dt <- mean(diff(DATA()$t))
        fnyq <- 1/dt/2
        kf <- req(input$kf) |> as.double()
        fmax <- round(kf*fnyq)
        
        k <- req(input$k) |> as.integer()
        .wl <- round(nrow(DATA())/k)
        .ovlp <- req(input$ovlp) |> as.double()
        .zp <- req(input$zp) |> as.double()
        .norm <- req(input$norm) |> as.logical()
        .osc <- req(input$osc) |> as.logical()
        .scale <- req(input$scale) |> as.logical()
        if(.norm==FALSE) .scale <- FALSE
        
        req(WAVE())
        PLOT <- seewave::spectro(
          WAVE(),
          flim=c(0,fmax/1000),
          wl=.wl,
          osc=.osc,
          ovlp=.ovlp,
          scale=.scale,
          zp=.zp,
          plot=TRUE,
          norm=.norm)
        PLOT
        
      })
    })
}


