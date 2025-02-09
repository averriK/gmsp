#' buildTS
#'
#' @param x data.table
#' @param dt numeric
#' @param ts numeric
#' @param Units character
#' @param Fmax integer
#' @param kNyq numeric
#' @param Resample boolean
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer
#' @param FlatZeros boolean
#' @param AstopAT numeric
#' @param ApassAT numeric
#' @param TrimZeros boolean
#' @param DetrendAT boolean
#' @param DetrendVT boolean
#' @param DetrendDT boolean
#' @param Normalize boolean
#' @param Output character
#'
#' @return list
#' @export 
#'
#' @examples
#'
#' @import data.table
#' @importFrom stats na.omit
#' @importFrom stats sd
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom seewave ffilter
#' @importFrom signal resample
#' @importFrom stringr str_split
#' @importFrom purrr map
#'
#'
build_TS <- function(
    x, ts=NULL ,dt=NULL, Units,
    DetrendAT=FALSE,
    DetrendVT=FALSE,
    DetrendDT=FALSE,
    Fmax = 16,
    kNyq=3.125, #>2.5
    Resample = TRUE,
    TargetUnits = "mm",
    NW = 128,
    OVLP = 75,
    FlatZeros=FALSE,
    AstopAT=1e-4,
    ApassAT=1e-3,
    TrimZeros=FALSE,
    Normalize = FALSE,
    Output=NULL) {
  on.exit(expr = {rm(list = ls())}, add = TRUE)
  . <- NULL
  X <- copy(x) |> as.data.table()
  
  NP <- nrow(X)
  NW_min <- NP/2
  if(NW>NW_min){
    L2 <- 2^floor(log2(NW_min))
    U2 <- 2^ceiling(log2(NW_min))
    if (abs(L2 - NW_min) < abs(U2 - NW_min)) {
      NW <- min(as.integer(L2),NP)
    } else {
      NW <- min(as.integer(U2),NP)
    }
  }
  if(!is.null(ts)){
    dts <- diff(ts)
    dt <- mean(dts)
  }
  
  # Check if dt is a rational number
  if(!.isRational(dt)){
    warning("Time step is not a rational number. Rounding to 3 decimal places (max sampling frequency 1kHz)")
    dt <- round(dt,3)
  }
  ts <- seq(0,(NP-1)*dt,by=dt)
  Fs <- 1 / dt
  
  ## Scale Units  ----
  SFU <- 1
  
  if (grepl(Units, pattern = "[///+]")) {
    Units <- (str_split(Units, pattern = "[///+]") |> unlist())[1]
  }
  if (!(tolower(Units) %in% c("mm", "cm", "m", "gal", "g"))) return(NULL)
  if (tolower(Units) != TargetUnits) {
    SFU <- .getSF(SourceUnits = tolower(Units), TargetUnits = TargetUnits)
    X <- X[,.(sapply(.SD, function(x) {x * SFU}))]
  }
  
  OCID <- names(X)
  
  # Export RAW record scaled to TargetUnits
  
  ATo <- data.table(ts=ts, Units=TargetUnits,X)
  if(!is.null(Output) && Output=="ATo"){return(ATo)}
  
  # setnames(RTSW,old=OCID,new=paste0("AT.",OCID))
  
  ## Scale record ----
  if(Normalize==TRUE){
    Ao <- apply(X, 2, function(x) { max(abs(x))})
    X <-X[, .(sapply(.SD, function(x){x/max(abs(x))}))]
  }
  
  # Detrend
  if(DetrendAT==TRUE){
    X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  
  
  ## Resample ----
  if(Resample){
    TargetFs <- as.integer(max(2.5,kNyq)*Fmax) # 80/200 Hz
    X <- .resample(X,dt=dt,TargetFs=TargetFs,Fmax=Fmax,NW=NW,OVLP=OVLP)
    Fs <- TargetFs
    dt <- 1/Fs
  }
  ## Case #2. Acceleration Time Histories ----
  if(FlatZeros==TRUE){
    Wo <- X[,.(sapply(.SD, function(x) {.taperA(x,Astop=SFU*AstopAT,Apass=SFU*ApassAT)}))]
  } else {
    Wo <- X[,.(sapply(.SD, function(x) {rep(1, length(x))}))]
  }
  
  AT <- X
  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  # AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  # Detrend
  if(DetrendAT==TRUE){
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  
  
  names(AT) <- OCID
  # browser()
  VT <- AT[, lapply(.SD, function(x){ .integrate(dx=x,dt=dt,NW=NW,OVLP=OVLP) })]
  
  if(nrow(VT) > nrow(Wo)){
    O <- data.table(sapply(Wo, function(x) rep(0, nrow(VT)-nrow(Wo))))
    Wo <- rbind(Wo, O)
  }
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  # VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  if(DetrendVT==TRUE){
    VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  
  
  names(VT) <- OCID
  
  DT <- VT[, lapply(.SD, function(x){ .integrate(dx=x,dt=dt,NW=NW,OVLP=OVLP) })]
  
  if(nrow(DT) > nrow(Wo)){
    O <- data.table(sapply(Wo, function(x) rep(0, nrow(DT)-nrow(Wo))))
    Wo <- rbind(Wo, O)
  }
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  
  ## Detrend DT ----
  if(DetrendDT==TRUE){
    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  names(DT) <- OCID
  
  
  
  
  ## Derivate DT ----
  VT <- DT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
  
  if(DetrendVT==TRUE){
    VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  ## Derivate VT ----
  AT <-  VT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
  
  ## Detrend AT ----
  if(DetrendAT==TRUE){
    AT <- AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  ## Taper Zeros ----
  # Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
  
  if(FlatZeros==TRUE){
    Wo <-AT[,.(sapply(.SD, function(x) {.taperA(x,Astop=SFU*AstopAT,Apass=SFU*ApassAT)}))]
  } else {
    Wo <- AT[,.(sapply(.SD, function(x) {rep(1, length(x))}))]
  }
  
  
  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  
  
  # Trim Zeros
  if(TrimZeros){
    idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
    AT <- AT[idx]
    VT <- VT[idx]
    DT <- DT[idx]
  }
  
  # Fix trend
  if(DetrendAT==TRUE){
    AT <- AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  if(DetrendVT==TRUE){
    VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  if(DetrendDT==TRUE){
    DT <- DT[, .(sapply(.SD, function(x){x-mean(x)}))]
  }
  
  
  
  ## Normalize Maximum ----
  if(Normalize==TRUE){
    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
    VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
    DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
  }
  
  names(AT) <- OCID
  names(VT) <- OCID
  names(DT) <- OCID
  if(!is.null(Output) && Output=="AT"){return(AT)}
  if(!is.null(Output) &&  Output=="VT"){return(VT)}
  if(!is.null(Output) && Output=="DT"){return(DT)}
  
  ## Summary ----
  NP <-  nrow(AT)
  ts <- seq(0,dt*(NP-1),dt)
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  ## Pack Time Series  ----
  
  TSW <- data.table(ts=ts, AT = AT, VT = VT, DT = DT)
  if(!is.null(Output) && Output=="TSW"){return(TSW)}
  
  ivars <- c("ts")
  mvars <- colnames(TSW[, -c("ts")])
  AUX <- data.table::melt(TSW, id.vars = ivars, measure.vars = mvars) |> na.omit()
  TSL <- AUX[,.(t=ts,s=value,ID=gsub("\\..*$", "", variable), OCID=gsub("^[^.]*\\.", "", variable))]
  if(!is.null(Output) && Output=="TSL"){return(TSL)}
  
  ## Trim Records,  ----
  # if(TrimZeros){
  #   #Warning: Different time scales from now on
  #   TSL <- TSL[,.trimZeros(.SD),by=c("OCID","ID")]}
  
  ## Return ----
  return(list(ATo=ATo,TSL=TSL,TSW=TSW,Wo=Wo,Ao=Ao,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = Units))
  
  
}
