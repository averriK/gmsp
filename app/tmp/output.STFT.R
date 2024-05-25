buildOutput.STFT <- function(id,ts){
  
  moduleServer(
    id=id, 
    module=function(input, output, session){
      
      kf <- req(input$kf) |> as.double()
      k <- req(input$k) |> as.integer()
      ovlp <- req(input$ovlp) |> as.double() 
      
      
      
      
      
      # STFT AT ----
      
      output$STFT.AT.H1 <- renderPlot({
        DT <- ts[ID=="AT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      output$STFT.AT.H2 <- renderPlot({
        DT <- ts[ID=="AT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      output$STFT.AT.UP <- renderPlot({
        DT <- ts[ID=="AT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      # STFT VT ----
      
      output$STFT.VT.H1 <- renderPlot({
        DT <- ts[ID=="VT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      output$STFT.VT.H2 <- renderPlot({
        DT <- ts[ID=="VT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      output$STFT.VT.UP <- renderPlot({
        DT <- ts[ID=="VT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      # STFT DT ----
      
      output$STFT.DT.H1 <- renderPlot({
        DT <- ts[ID=="DT"& OCID=="H1"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      output$STFT.DT.H2 <- renderPlot({
        DT <- ts[ID=="DT"& OCID=="H2"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
      output$STFT.DT.UP <- renderPlot({
        DT <- ts[ID=="DT"& OCID=="UP"]
        dt <- mean(diff(DT$t))
        fnyq <- 1/dt/2
        fmax <- round(kf*fnyq)
        n <- nrow(DT)
        wl <- round(n/k)
        STFT <- .getSTFT(s=DT$s,dt=dt)
        seewave::spectro(STFT,flim=c(0,fmax/1000),wl=wl,osc=TRUE,dB=NULL,ovlp=ovlp,zp=256,plot=TRUE,norm=TRUE)
      })
      
    })
}
