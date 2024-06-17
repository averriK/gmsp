
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


SDOF.server <- function(id,.data,color.palette="Dynamic"){
  moduleServer(
    id, 
    function(input, output, session){
      
      DATA <- reactive({
        req(input$ocid)
        req(.data)
        stopifnot(length(input$ocid)==1&& input$ocid %in% c("H1","H2","UP"))
        DT <- .data[ID =="AT" & OCID == input$ocid,.(t,s)]
        if(nrow(DT)==0) return(NULL)
        xi <- req(input$xi)
        DT <- gmsp::build_PSA(.x=DT,Units=UNITS,xi=xi) #UNITS is a global variable
        
        
        req(input$smoothing)
        window <- req(input$window) |> as.double()
        cv <- req(input$cv) |> as.logical()
        sigma <- req(input$sigma) |> as.double()
        
        DATA <- data.table()
        AUX <- gmsp::smooth_Spectra(fs=DT$Tn,A=DT$PSA,method=input$smoothing,window=window,cv=cv,sigma=sigma)
        DATA <- rbindlist(list(DATA, data.table(ID="PSA",X=AUX$fs,Y=AUX$As)))
        
        AUX <- gmsp::smooth_Spectra(fs=DT$Tn,A=DT$PSV,method=input$smoothing,window=window,cv=cv,sigma=sigma)
        DATA <-  rbindlist(list(DATA, data.table(ID="PSV",X=AUX$fs,Y=AUX$As)))
        
        AUX <- gmsp::smooth_Spectra(fs=DT$Tn,A=DT$SD,method=input$smoothing,window=window,cv=cv,sigma=sigma)
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


STFT.server <- function(id,.data,series){
  moduleServer(
    id, 
    function(input, output, session){
      stopifnot(series %in% c("AT","VT","DT"))
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


CWT.server <- function(id,.data,series){
  moduleServer(
    id, 
    function(input, output, session){
      stopifnot(series %in% c("AT","VT","DT"))
    
      
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


