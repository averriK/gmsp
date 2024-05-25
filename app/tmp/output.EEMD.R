buildOutput.EEMD <- function(id,ts){
  
  hht::EEMD(sig=DT$s, tt=DT$t,nimf=nimf,max.imf=max.imf,boundary=boundary,noise.amp=noise.amp, noise.type=noise.type,trials=trials,stop.rule=stop.rule,trials.dir = DIR,verbose=FALSE)
  AUX <- EEMDCompile(trials.dir = DIR, trials=trials, nimf=nimf) |> suppressWarnings()
  
  
  
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
      
      output$EEMD.AT.H1 <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      output$EEMD.AT.H2 <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      output$EEMD.AT.UP <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      # VT ----
      
      output$EEMD.VT.H1 <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      output$EEMD.VT.H2 <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      output$EEMD.VT.UP <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      # DT ----
      
      output$EEMD.DT.H1 <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      output$EEMD.DT.H2 <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
      
      output$EEMD.DT.UP <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,dt=dt,method="eemd",boundary=boundary, max.imf=max.imf,noise.type=noise.type,noise.amp=noise.amp,trials=trials,stop.rule=stop.rule,plot=TRUE)
        DATA <- AUX$plot.data
        nimf <- AUX$nimf
        plot.highchart(
          color.palette ="ag_Sunset",
          plot.height = max(500,150*nimf),
          plot.type="line",
          legend.layout="horizontal",
          legend.show=FALSE,
          xAxis.log=FALSE,yAxis.log=FALSE,
          yAxis.legend="IMF(t)",xAxis.legend="t",
          data=DATA)
      })
    })
}
