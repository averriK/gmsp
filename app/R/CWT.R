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
  
  
} 

CWT.mainPanel <- buildMainPanel(  id="CWT",render="plotOutput",
                                  title="Continuous Wavelet Transform (CWT)",
                                  title.AT="Acceleration (AT)",
                                  title.VT="Velocity (VT)",
                                  title.DT="Displacement (DT)",
                                  height="700px")

CWT.tabPanel <- buildTabPanel(id="CWT")




CWT.server <- function(id,.data,series){
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
        stopifnot(length(input$ocid)==1&& input$ocid %in% c("H1","H2","UP"))
        req(.data)
        DT <- .data[ID %in% series & OCID %in% input$ocid,.(t,s=as.matrix(s))]
        
        if(nrow(DT)==0) return(NULL)
        dt <- mean(diff(DT$t))
        mother <- req(input$mother)
        j <- req(input$j) |> as.integer()
        dj <- 1/j
        DT <- biwavelet::wt(DT, mother = mother,dt=dt,dj=dj,do.sig=FALSE)
        DT
      })
      
      output[[series]] <- renderPlot({
        req(DATA())
        phase <- req(input$phase)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        PLOT <- plot(DATA(), plot.cb = TRUE, plot.phase = phase)
        PLOT
        
      })
    })
}


