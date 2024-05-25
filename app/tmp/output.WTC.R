output.WTC <- function(id,ts){
  
  
  moduleServer(
    id=id, 
    module=function(input, output, session){
      
      mother <- input$mother
      j <- req(input$j) |> as.integer()
      dj <- 1/j
      phase <- req(input$phase) 
      if(is.null(phase)){phase <- FALSE}
      
      # WTC AT ----
      
      output$WTC.AT.H1 <-  renderPlot({
        DT <- ts[ID=="AT"& OCID=="H1"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      output$WTC.AT.H2 <-  renderPlot({
        DT <- ts[ID=="AT"& OCID=="H2"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      output$WTC.AT.UP <-  renderPlot({
        DT <- ts[ID=="AT"& OCID=="UP"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      # WTC VT ----
      
      output$WTC.VT.H1 <-  renderPlot({
        DT <- ts[ID=="VT"& OCID=="H1"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      output$WTC.VT.H2 <-  renderPlot({
        DT <- ts[ID=="VT"& OCID=="H2"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      output$WTC.VT.UP <-  renderPlot({
        DT <- ts[ID=="VT"& OCID=="UP"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      # WTC DT ----
      
      output$WTC.DT.H1 <-  renderPlot({
        DT <- ts[ID=="DT"& OCID=="H1"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      output$WTC.DT.H2 <-  renderPlot({
        DT <- ts[ID=="DT"& OCID=="H2"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      output$WTC.DT.UP <-  renderPlot({
        DT <- ts[ID=="DT"& OCID=="UP"]
        WTC <- .getWTC(s=DT$s,t=DT$t,dj=dj,mother=mother)
        par(oma = c(0, 0, 0, 1), mar = c(5, 4, 4, 5) + 0.2)
        plot(WTC, plot.cb = TRUE, plot.phase = phase)
      })
      
      
      
    })
}
