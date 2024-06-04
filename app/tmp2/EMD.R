EMD.sidebar <- function(id){
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


EMD.server <- function(id,.data,series,yAxis.legend="s(t)",xAxis.legend="t",color.palette="Dynamic"){
  moduleServer(
    id, 
    function(input, output, session){
      stopifnot(series %in% c("AT","VT","DT"))
      observe({
        cat("Current OCID value:", input$ocid, "\n")  
        cat("Current max.imf value:", input$max.imf, "\n")  
        cat("Current boundary value:", input$boundary, "\n")  
        cat("Current stop.rule value:", input$stop.rule, "\n")  
      })
      
      DATA <- reactive({
        req(input$ocid)
        stopifnot(length(input$ocid)==1&& input$ocid %in% c("H1","H2","UP"))
        req(.data)
        DT <- .data[ID %in% series & OCID %in% input$ocid,.(t,s)]
        if(nrow(DT)==0) return(NULL)
        dt <- mean(diff(DT$t))
        req(input$max.imf)
        req(input$stop.rule)
        req(input$boundary)
        DT <- gmsp::buildIMF(DT$s,dt=dt,method="emd",boundary=input$boundary, max.imf=input$max.imf,stop.rule=input$stop.rule,plot=TRUE,verbose=TRUE)
        DT
      })
      
      output[[series]] <- renderHighchart({
        req(DATA())
        IMF <- DATA()$plot.data
        PLOT <- buildPlot::buildPlot(
          library="hc",
          color.palette = color.palette,
          plot.type="line",
          legend.layout="horizontal",
          legend.show=TRUE,
          yAxis.label=FALSE,
          group.legend="IMF",
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend=yAxis.legend,
          xAxis.legend=xAxis.legend,
          data=IMF)
        PLOT
        
      })
    })
}


