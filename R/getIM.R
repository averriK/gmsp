
#' Title
#'
#' @param a data.table Time Series of acceleration records
#' @param v data.table Time Series of velocity records
#' @param d data.table Time Series of displacement records
#' @param dt numeric Time step
#' @param UN character Units of source records
#' @param TargetUnits character Units of target records
#'
#' @return data.table
#' @export getIM
#'
#' @examples
#'
#' @import data.table
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom signal resample
#' @importFrom pracma detrend
#' @importFrom stringr str_split
#' @importFrom purrr map

getIM <- function(a=NULL,v=NULL,d=NULL,dt=NULL,UN=NULL,TargetUnits="mm"){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  # Check internal  -------------------------------------------------------------
  OK <- !is.null(a) && !is.null(dt) && !is.null(UN)
  stopifnot(OK)
  AT <- copy(a)
  VT <- copy(v)
  DT <- copy(d)

  # Check Records -------------------------------------------------------------

  if(
    nrow(AT)==0||
    dt==0||
    any(is.na(AT))||
    max(abs(AT))==0 #||
    #any(sapply(AT, function(x){max(abs(x))})==0)
  ){
    # Null Record
    return(NULL)
  }
  IM <- data.table()

  OCID <- colnames(AT)

  if(length(unique(OCID))<3){
    # Invalid Header
    return(NULL)
  }

  IM[,OCID:=OCID]

  if(grepl(UN,pattern = "[///+]")){
    UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
  }
  if(!(tolower(UN) %in% c("mm","cm","m","gal","g"))){
    # Invalid Units
    return(NULL)
  }

  NP <- nrow(AT)
  tn <- seq(from=0,length.out=NP,by=dt)
  Fs=1/dt
  # Scale Records -------------------------------------------------------------
  g <- .getG(TargetUnits) #GMSP$g
  if(UN!=TargetUnits) {
    SFU <- .getSF(UN,TargetUnits = TargetUnits)
    AT[, (OCID) := lapply(.SD, function(x){SFU*x}), .SDcols=OCID]
    IM[,UN:=TargetUnits]
  } else {
    SFU <- 1
    IM[,UN:=UN]
  }
  IM[,SFU:=SFU]
  IM[,Scaled:=FALSE]
  # fs  <- Fs/2*.linspace(0,1,NUP) # ?
  IM[,NP:=NP]
  IM[,dt:=dt]
  IM[, Fs:= Fs]


  # AT Intensity  -------------------------------------------------------------
  g <- .getG(TargetUnits) #GMSP$g

  IM[,IA := AT[, sapply(.SD, function(x){.getIA(x,dt=dt,g=g)})]] # Arias Intensity
  IM[,IAu := AT[, sapply(.SD, function(x){.getIA(max(x,0),dt=dt,g=g)})]] # Arias Intensity
  IM[,IAd := AT[, sapply(.SD, function(x){.getIA(min(x,0),dt=dt,g=g)})]] # Arias Intensity


  IM[,PGA :=AT[, sapply(.SD, function(x){max(abs(x))})]] #  Peak values
  IM[,AZC := AT[, sapply(.SD,   .getZC)]]
  IM[,ARMS := AT[, sapply(.SD,  .getRMS)]]
  IM[,CAV5 := AT[, sapply(.SD,  function(x){.getCAV5(x,tn=tn,UN=TargetUnits)})]]
  IM[,TmA :=  AT[, sapply(.SD, function(x){.getTm(X=x,Fs=Fs,fmin=0.25,fmax=15)})]]
  IM[,CRC32 := AT[, sapply(.SD, .getCRC)]]
  IM[,ATo :=AT[, sapply(.SD, function(x){first(x)})]] #  Initial values
  IM[,ATn :=AT[, sapply(.SD, function(x){last(x)})]] #  Last values
  # VT Intensity  -------------------------------------------------------------
  if(!is.null(VT)){

    IM[,PGV :=VT[, sapply(.SD, function(x){max(abs(x))})]] #  Peak values
    IM[,VZC := VT[, sapply(.SD,   .getZC)]]
    IM[,VRMS := VT[, sapply(.SD,  .getRMS)]]
    IM[,VTo := VT[, sapply(.SD, function(x){first(x)})]] #  Initial values
    IM[,VTn := VT[, sapply(.SD, function(x){last(x)})]] #  Last values
  }
  else {
    IM[,PGV :=NA] #  Peak values
    IM[,VZC := NA]
    IM[,VRMS := NA]
    IM[,VTo := NA] #  Initial values
    IM[,VTn := NA] #  Last values
  }

  # DT Intensity  -------------------------------------------------------------
  if(!is.null(DT)){
    IM[,PGD :=DT[, sapply(.SD, function(x){max(abs(x))})]] #  Peak values
    IM[,DZC := VT[, sapply(.SD,   .getZC)]]
    IM[,DRMS := DT[, sapply(.SD,  .getRMS)]]
    IM[,DTo := DT[, sapply(.SD, function(x){first(x)})]] #  Initial values
    IM[,DTn := DT[, sapply(.SD, function(x){last(x)})]] #  Last values
  } else {
    IM[,PGD :=NA] #  Peak values
    IM[,DZC := NA]
    IM[,DRMS := NA]
    IM[,DTo := NA] #  Initial values
    IM[,DTn := NA] #  Last values
  }

  # Newmark Displacements  -------------------------------------------------------------
  #
  IM[,DN05 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.05,g=g)})]]
  IM[,DN10 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.10,g=g)})]]
  IM[,DN15 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.15,g=g)})]]
  IM[,DN20 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.20,g=g)})]]
  IM[,DN25 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.25,g=g)})]]
  IM[,DN30 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.30,g=g)})]]
  IM[,DN40 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.40,g=g)})]]
  IM[,DN50 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.50,g=g)})]]
  IM[,DN75 := AT[, sapply(.SD, function(x){.getDN(x,tn=tn,kh=0.75,g=g)})]]

  # Bracketed (Hussid) Duration
  IM[, Dmax := tn |> last()]
  IM[, D0595 := AT[, sapply(.SD, function(x){.getHBD(x,tn,a=0.05,b=0.95,g=g)})]]
  IM[, D2080 := AT[, sapply(.SD, function(x){.getHBD(x,tn,a=0.20,b=0.80,g=g)})]]
  IM[, D0575 := AT[, sapply(.SD, function(x){.getHBD(x,tn,a=0.05,b=0.75,g=g)})]]

  # Damage Indices
  IM[,PPI:=sqrt((ARMS^3)*D0595)]
  IM[,EPI:= 0.9/pi*IA*2*g*D0595]
  IM[,PDI:= ifelse(AZC>0,IA*(Dmax/AZC)^2,0)]

  return(IM)
}




