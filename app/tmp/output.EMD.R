output.EMD <- function(id,ts){
  
  # AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf,stoprule=stop.rule)
  
  
  
  moduleServer(
    id=id, 
    module=function(input, output, session){
      
      boundary <- req(input$boundary)
      max.imf <- req(input$max.imf) |> as.integer()
      stop.rule <- req(input$stop.rule)
      # AT ----
      
      output$EMD.AT.H1 <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE,verbose=TRUE)
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
      
      output$EMD.AT.H2 <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.AT.UP <- renderHighchart({
        DT <- ts[ID=="AT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.VT.H1 <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.VT.H2 <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.VT.UP <- renderHighchart({
        DT <- ts[ID=="VT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.DT.H1 <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.DT.H2 <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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
      
      output$EMD.DT.UP <- renderHighchart({
        DT <- ts[ID=="ADT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        AUX <- buildIMF(s=DT$s,t=DT$t,method="emd",boundary=boundary, max.imf=max.imf,stop.rule=stop.rule,plot=TRUE)
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