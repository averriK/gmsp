TS.server <- function(id,.data,series,yAxis.legend="s(t)",xAxis.legend="t",color.palette="Dynamic"){
  moduleServer(
    id, 
    function(input, output, session){
      observe({
        cat("Current OCID value:", input$ocid, "\n")  
      })
      stopifnot(series %in% c("AT","VT","DT"))
      DATA <- reactive({
        req(input$ocid)
        req(.data)
        stopifnot(length(input$ocid)==1&& input$ocid %in% c("H1","H2","UP"))
        DT <- .data[ID %in% series & OCID %in% input$ocid,.(ID=OCID,X=t,Y=s)]
        if(nrow(DT)==0) return(NULL)
        DT
      })
      output[[series]] <- renderHighchart({
        req(DATA())
        PLOT <-  buildPlot::buildPlot(
          library="hc",
          color.palette = color.palette,
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend=yAxis.legend,
          xAxis.legend=xAxis.legend,
          data=DATA())
        PLOT
        
      })
    })
}

TS.sidebar <- function(id) {
  ns <- NS(id)
  tagList(
    wellPanel(
      titlePanel("OCID"),
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
  )
}



