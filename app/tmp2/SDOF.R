SDOF.navbarPageUI <- function(id,title="Single Degree of Freedom (SDOF)",title.AT="Pseudo Spectral Acceleration (PSA)",title.VT="Pseudo Spectral Velocity (PSV)",title.DT="Spectrl Displacement (SD)",height="500px"){
  ns <- NS(id)
  navbarPage(title=title,
             id = id,
             fluidRow(
               column(
                 width=8,
                 tabsetPanel(
                   tabPanel(title=title.AT,highchartOutput(ns("AT"),height = height)),
                   tabPanel(title=title.VT,highchartOutput(ns("VT"),height = height)),
                   tabPanel(title=title.DT,highchartOutput(ns("DT"),height = height))
                 )
               )
             )
  ) # navbarPage
}


SDOF.sidebar <- function(id){
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
    ),
    
    
    
    
    
    
    
  )
}

SDOF.server <- function(id,.data,color.palette="Dynamic"){
  moduleServer(
    id, 
    function(input, output, session){
      observe({
        cat("Current xi value:", input$xi, "\n")
        cat("Current OCID value:", input$ocid, "\n")  
        cat("Current smoothing value:", input$smoothing, "\n")  
        cat("Current cv value:", input$cv, "\n")  
        cat("Current sigma value:", input$sigma, "\n")  
        cat("Current window value:", input$window, "\n")  
      })
      
      DATA <- reactive({
        req(input$ocid)
        req(.data)
        stopifnot(length(input$ocid)==1&& input$ocid %in% c("H1","H2","UP"))
        DT <- .data[ID =="AT" & OCID == input$ocid,.(t,s)]
        if(nrow(DT)==0) return(NULL)
        xi <- req(input$xi)
        DT <- gmsp::buildPSA(.x=DT,Units=UNITS,xi=xi) #UNITS is a global variable
        
        
        req(input$smoothing)
        window <- req(input$window) |> as.double()
        cv <- req(input$cv) |> as.logical()
        sigma <- req(input$sigma) |> as.double()
        
        DATA <- data.table()
        AUX <- gmsp::smoothSpectra(fs=DT$Tn,A=DT$PSA,method=input$smoothing,window=window,cv=cv,sigma=sigma)
        DATA <- rbindlist(list(DATA, data.table(ID="PSA",X=AUX$fs,Y=AUX$As)))
        
        AUX <- gmsp::smoothSpectra(fs=DT$Tn,A=DT$PSV,method=input$smoothing,window=window,cv=cv,sigma=sigma)
        DATA <-  rbindlist(list(DATA, data.table(ID="PSV",X=AUX$fs,Y=AUX$As)))
        
        AUX <- gmsp::smoothSpectra(fs=DT$Tn,A=DT$SD,method=input$smoothing,window=window,cv=cv,sigma=sigma)
        DATA <-  rbindlist(list(DATA, data.table(ID="SD",X=AUX$fs,Y=AUX$As)))
        
        DATA
      })
      
      output[["AT"]] <- renderHighchart({
        req(DATA())
        PLOT <-  buildPlot::buildPlot(
          library="hc",
          color.palette = color.palette,
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=TRUE,yAxis.log=FALSE,
          xAxis.legend="Tn[s]",
          yAxis.legend=paste0("PSA[",UNITS,"/s2]"),
          data=DATA()[ID=="PSA"])
        PLOT
      })
      
      
      
      output[["VT"]] <- renderHighchart({
        req(DATA())
        PLOT <-  buildPlot::buildPlot(
          library="hc",
          color.palette = color.palette,
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=TRUE,yAxis.log=FALSE,
          xAxis.legend="Tn[s]",
          yAxis.legend=paste0("PSV[",UNITS,"/s]"),
          data=DATA()[ID=="PSV"])
        PLOT
        
      })
      
      output[["DT"]] <- renderHighchart({
        req(DATA())
        PLOT <-  buildPlot::buildPlot(
          library="hc",
          color.palette = color.palette,
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=TRUE,yAxis.log=FALSE,
          xAxis.legend="Tn[s]",
          yAxis.legend=paste0("SD[",UNITS,"]"),
          data=DATA()[ID=="SD"])
        PLOT
        
      })
      
      
    })
}



