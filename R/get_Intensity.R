#' Title
#'
#' @param TSL data.table
#' @param TargetUnits character
#'
#' @return
#' @export
#'
#' @examples



get_Intensity <- function(TSL,TargetUnits="mm"){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  . <- NULL
  # Arias Intensity -------------------------------------------------------------
  g <- .getG(TargetUnits) #GMSP$g
  IA <- TSL[ID=="AT",.(ID="AI",value=.getAI(x=s,t=t,g=g)),by=.(RecordSN,DIR,OCID)]
  IAu <- TSL[ID=="AT",.(ID="AIu",value=.getAI(x=max(s,0),t=t,g=g)),by=.(RecordSN,DIR,OCID)]
  IAd <- TSL[ID=="AT",.(ID="AId",value=.getAI(x=min(s,0),t=t,g=g)),by=.(RecordSN,DIR,OCID)]
  
  # Peak Values  -------------------------------------------------------------
  PGA <- TSL[ID=="AT",.(ID="PGA",value=.getPeak(s)),by=.(RecordSN,DIR,OCID)]
  PGV <- TSL[ID=="VT",.(ID="PGV",value=.getPeak(s)),by=.(RecordSN,DIR,OCID)]
  PGD <- TSL[ID=="DT",.(ID="PGD",value=.getPeak(s)),by=.(RecordSN,DIR,OCID)]
  
  # RMS Values  -------------------------------------------------------------
  ARMS <- TSL[ID=="AT",.(ID="ARMS",value=.getRMS(s)),by=.(RecordSN,DIR,OCID)]
  VRMS <- TSL[ID=="VT",.(ID="VRMS",value=.getRMS(s)),by=.(RecordSN,DIR,OCID)]
  DRMS <- TSL[ID=="DT",.(ID="DRMS",value=.getRMS(s)),by=.(RecordSN,DIR,OCID)]
  
  
  # Newmark Displacements  -------------------------------------------------------------
  kh_values <- c(0.01,0.02,0.05, 0.10, 0.15, 0.20, 0.25, 0.30,0.35, 0.40,0.45, 0.50,0.55,0.60,0.65,0.70,0.75)
  NDT <- rbindlist(
    sapply(kh_values, function(kh) {
      TSL[ID == "AT", .(
        ID = sprintf("DN%02d", round(kh * 100)), 
        value = gmsp::get_Newmark(AT = s, t = t, kh = kh, FULL = FALSE)), 
        by = .(RecordSN, DIR, OCID)]}, simplify = FALSE))
  
  # Sample properties -------------------------------------------------------------
  NP <- TSL[ID=="AT",.(ID="NP",value=length(s)),by=.(RecordSN,DIR,OCID)]
  dt <- TSL[ID=="AT",.(ID="dt",value=mean(diff(t))),by=.(RecordSN,DIR,OCID)]
  Fs <- TSL[ID=="AT",.(ID="Fs",value=1/mean(diff(t))),by=.(RecordSN,DIR,OCID)]
  
  # Initial Values -------------------------------------------------------------
  ATo <- TSL[ID=="AT",.(ID="ATo",value=first(s)),by=.(RecordSN,DIR,OCID)]
  VTo <- TSL[ID=="VT",.(ID="VTo",value=first(s)),by=.(RecordSN,DIR,OCID)]
  DTo <- TSL[ID=="DT",.(ID="DTo",value=first(s)),by=.(RecordSN,DIR,OCID)]
  
  # Last Values -------------------------------------------------------------
  ATn <- TSL[ID=="AT",.(ID="ATn",value=last(s)),by=.(RecordSN,DIR,OCID)]
  VTn <- TSL[ID=="VT",.(ID="VTn",value=last(s)),by=.(RecordSN,DIR,OCID)]
  DTn <- TSL[ID=="DT",.(ID="DTn",value=last(s)),by=.(RecordSN,DIR,OCID)]
  # Zero Crossing -------------------------------------------------------------
  AZC <- TSL[ID=="AT",.(ID="AZC",value=.getZC(s)),by=.(RecordSN,DIR,OCID)]
  VZC <- TSL[ID=="VT",.(ID="VZC",value=.getZC(s)),by=.(RecordSN,DIR,OCID)]
  DZC <- TSL[ID=="DT",.(ID="DZC",value=.getZC(s)),by=.(RecordSN,DIR,OCID)]
  # Bracketed (Hussid) Duration ----------------------------------------------
  Dmax <- TSL[ID=="AT",.(ID="Dmax",value=last(t)),by=.(RecordSN,DIR,OCID)]
  D0595 <- TSL[ID=="AT",.(ID="D0595",value=.getHBD(s,t=t,a=0.05,b=0.95,g=g)),by=.(RecordSN,DIR,OCID)]
  D2080 <- TSL[ID=="AT",.(ID="D2080",value=.getHBD(s,t=t,a=0.20,b=0.80,g=g)),by=.(RecordSN,DIR,OCID)]
  D0575 <- TSL[ID=="AT",.(ID="D0575",value=.getHBD(s,t=t,a=0.05,b=0.75,g=g)),by=.(RecordSN,DIR,OCID)]
  
  # Tm -------------------------------------------------------------
  TmA <- TSL[ID=="AT",.(ID="TmA",value=.getTm(s,t=t,fmin=0.1,fmax=25)),by=.(RecordSN,DIR,OCID)]
  TmV <- TSL[ID=="VT",.(ID="TmV",value=.getTm(s,t=t,fmin=0.1,fmax=25)),by=.(RecordSN,DIR,OCID)]
  TmD <- TSL[ID=="DT",.(ID="TmD",value=.getTm(s,t=t,fmin=0.1,fmax=25)),by=.(RecordSN,DIR,OCID)]
  
  # Build Intensity Table
  ITL <- rbindlist(list(IA,IAu,IAd,PGA,PGV,PGD,ARMS,VRMS,DRMS,NDT,NP,dt,Fs,ATo,VTo,DTo,ATn,VTn,DTn,AZC,VZC,DZC,Dmax,D0595,D2080,D0575,TmA,TmV,TmD), use.names = TRUE)
  ITW <- dcast(ITL, RecordSN + OCID + DIR ~ ID, value.var = "value")
  return(list(ITL=ITL,ITW=ITW))
}


.getG <- function(TargetUnits,g= 9806.650) {
  switch (tolower(TargetUnits),"m" = g/1000,"cm" = g/10,"mm" = g/1,NULL)
}

.getAI <- function(x,t,g=NULL,TargetUnits="mm"){
  if(is.null(g) & !is.null(TargetUnits)){g <- .getG(TargetUnits)}
  dt <- mean(diff(t))
  as.numeric(x %*% x)*dt*pi/(2*g)
}

.getRMS <- function(x){
  sqrt(1/length(x)*as.numeric(x%*%x))}

.getPeak <- function(x){
  max(abs(x))
}


.getZC <- function(x){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  NP <- length(x)
  ZC <- 0
  if(NP>2){ZC <- sum(sign( x[2:NP])== -sign(x[1:(NP-1)]))}
  return(ZC)
}

.getCRC <- function(x){toupper(digest(object=x,algo="crc32"))}

.getHBD <- function(x,t,a=0.05,b=0.95,g=9806.650) {
  on.exit(expr={rm(list = ls())}, add = TRUE)
  dt <-mean(diff(t))
  IA <- dt*(x%*%x)*pi/(2*g)
  SumIA <- pi*dt*cumsum(x^2)
  A <- a*2*g*IA
  B <- b*2*g*IA
  k <- 1
  while (SumIA[k]<A){
    k <- k+1
  }
  ta <- t[k]
  while (SumIA[k]<B){
    k <- k+1
  }
  tb <- t[k]
  D  <- tb-ta
  return(D)
}
.getTm <- function(x,t,fmin = 0,fmax = Inf){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  if(max(abs(x))==0) {return(0)}
  NP <- length(x)
  dt <- mean(diff(t))
  Fs <- 1/dt
  # NFFT <- nextn(NP,factors = 2)
  AW <- 1/NP*fft(x ,inverse = FALSE)
  df <- Fs/NP
  NUP=ceiling(NP/2)+1
  fs <- seq(from=1,to=NUP)*df
  fmin <- max(fmin,min(fs))
  fmax <- min(fmax,max(fs))
  # idx <- inrange(fs[1:NUP],lower=max(fmin,min(fs)),upper= min(fmax,max(fs)))
  idx <- fs>=fmin & fs<=fmax
  Co   <- sqrt(Re(AW[1:NUP]*Conj(AW[1:NUP])))
  # f1 <- max(fmin,min(fs))
  # f2 <- min(fmax,max(fs))
  # ix <- fs[1:NUP] %inrange% c(f1,(f2+fs[2]))
  Tm <-sum(Co[idx]^2/fs[idx])/sum(Co[idx]^2)
  return(Tm)
}

