

FFT.sidebar <- function(id){
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


FFT.server <- function(id,.data,series,yAxis.legend="A(fs)",xAxis.legend="fs",color.palette="Dynamic"){
  moduleServer(
    id, 
    function(input, output, session){
      observe({
        cat("Current OCID value:", input$ocid, "\n")  
        cat("Current smoothing value:", input$smoothing, "\n")  
        cat("Current kf value:", input$kf, "\n")  
        cat("Current zp value:", input$zp, "\n")  
        cat("Current cv value:", input$cv, "\n")  
        cat("Current sigma value:", input$sigma, "\n")  
        cat("Current window value:", input$window, "\n")  
      })
      stopifnot(series %in% c("AT","VT","DT"))
      DATA <- reactive({
        req(input$ocid)
        req(.data)
        stopifnot(length(input$ocid)==1&& input$ocid %in% c("H1","H2","UP"))
        DT <- .data[ID %in% series & OCID %in% input$ocid]
        if(nrow(DT)==0) return(NULL)
        zp <- req(input$zp) |> as.double()
        kf <- req(input$kf) |> as.double()
        DT <- gmsp::buildFFT(s=DT$s,t=DT$t,kf=kf,zp=zp)
        
        req(input$smoothing)
        window <- req(input$window) |> as.double()
        cv <- req(input$cv) |> as.logical()
        sigma <- req(input$sigma) |> as.double()
        AUX <- gmsp::smoothSpectra(fs=DT$fs,A=DT$A,method=input$smoothing,window=window,cv=cv,sigma=sigma)
        DATA <- data.table(ID=input$ocid,X=AUX$fs,Y=AUX$As)
        DATA
      })
      output[[series]] <- renderHighchart({
        req(DATA())
        PLOT <- buildPlot::buildPlot(
          library="hc",
          color.palette = color.palette,
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=TRUE,yAxis.log=FALSE,
          yAxis.legend=yAxis.legend,
          xAxis.legend=xAxis.legend,
          data=DATA())
        PLOT
        
      })
    })
}



