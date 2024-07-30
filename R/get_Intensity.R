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
  
  # Damage Indices -------------------------------------------------------------
  AUX <- ARMS[,.(RecordSN,DIR,OCID,ARMS=value)][D0595[,.(RecordSN,DIR,OCID,D0595=value)],on=.(RecordSN,DIR,OCID)]
  PPI <- AUX[,.(RecordSN,DIR,OCID,ID="PPI",value=sqrt((ARMS^3)*D0595))]
  
  AUX <- IA[,.(RecordSN,DIR,OCID,IA=value)][D0595[,.(RecordSN,DIR,OCID,D0595=value)],on=.(RecordSN,DIR,OCID)]
  EPI <- AUX[,.(RecordSN,DIR,OCID,ID="EPI",value=0.9/pi*IA*2*g*D0595)]
  
  
  AUX <- AZC[,.(RecordSN,DIR,OCID,AZC=value)][Dmax[,.(RecordSN,DIR,OCID,Dmax=value)],on=.(RecordSN,DIR,OCID)][IA[,.(RecordSN,DIR,OCID,IA=value)],on=.(RecordSN,DIR,OCID)]
  PDI <- AUX[,.(RecordSN,DIR,OCID,ID="PDI",value=IA*(Dmax/AZC)^2)]
  
  
  
  # Build Intensity Table
  ITL <- rbindlist(list(IA,IAu,IAd,PGA,PGV,PGD,ARMS,VRMS,DRMS,NP,dt,Fs,ATo,VTo,DTo,ATn,VTn,DTn,AZC,VZC,DZC,Dmax,D0595,D2080,D0575,TmA,TmV,TmD,PPI,EPI,PDI), use.names = TRUE)
  ITW <- dcast(ITL, RecordSN + OCID + DIR ~ ID, value.var = "value")
  return(list(ITL=ITL,ITW=ITW))
}





