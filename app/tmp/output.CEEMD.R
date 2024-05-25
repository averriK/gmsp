buildOutput.CEEMD <- function(id,ts){
  # hht::CEEMD(sig=DT$s, tt=DT$t,noise.amp=noise.amp, noise.type=noise.type,trials=trials,stop.rule=stop.rule,verbose=FALSE)
  
  
  moduleServer(
    id=id, 
    module=function(input, output, session){
      boundary <- req(input$boundary)
      max.imf <- req(input$max.imf) |> as.integer()
      stop.rule <- req(input$stop.rule)
      noise.type <- req(input$noise.type)
      noise.amp <- req(input$noise.amp) |> as.double()
      trials <- req(input$trials) |> as.integer()
      # AT ----
      
      output$CEEMD.AT.H1 <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s,dt=dt,method="ceemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE,verbose=TRUE)
        
      })
      
      # VT ----
      
      output$CEEMD.VT.H1 <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
      })
      
      # DT ----
      
      output$CEEMD.DT.H1 <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
      })
      
    })
}
