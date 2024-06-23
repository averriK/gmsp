#' buildTS
#'
#' @param x data.table
#' @param dt numeric
#' @param ts numeric
#' @param OrderTS integer
#' @param Units character
#' @param Fmax integer
#' @param kNyq numeric
#' @param Resample boolean
#' @param LowPass boolean
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer
#' @param Astop.AT numeric
#' @param Apass.AT numeric
#' @param Scale character

#'
#' @return list
#' @export 
#'
#' @examples
#'
#' @import data.table
#' @importFrom stats na.omit
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
    OrderTS=2, #0 Displacement, 1 Velocity, 2 Acceleration
    Fmax = 16,
    kNyq=3.125, #>2.5
    Resample = TRUE,
    LowPass = TRUE,
    TargetUnits = "mm",
    NW = 1024,
    OVLP = 75,
    Astop.AT=1e-4,
    Apass.AT=1e-3,
    Scale = "relative") {
  on.exit(expr = {rm(list = ls())}, add = TRUE)
  . <- NULL
  X <- copy(x) |> as.data.table()

  NP <- nrow(X)
  stopifnot(NP>4 && NP>=NW)
  if(!is.null(ts)){
    dt <- diff(ts) |> mean()
  }
  ts <- seq(0,(NP-1)*dt,by=dt)
  Fs <- 1 / dt
  
  

  ## Scale Units  ----

  if (grepl(Units, pattern = "[///+]")) {
    Units <- (str_split(Units, pattern = "[///+]") |> unlist())[1]
  }
  if (!(tolower(Units) %in% c("mm", "cm", "m", "gal", "g"))) return(NULL)

  if (tolower(Units) != TargetUnits) {
    SFU <- .getSF(SourceUnits = tolower(Units), TargetUnits = TargetUnits)
    # AT <- map(AT,function(x){x*SFU})
    # ATo[, (colnames(ATo)) := lapply(.SD, function(x) {x * SFU})]
    X <- X[,.(sapply(.SD, function(x) {x * SFU}))]
  }

  OCID <- names(X)
  
  # Export RAW record scaled to TargetUnits
  
  RTSW <- data.table(ts=ts, Units=TargetUnits,X)
  setnames(RTSW,old=OCID,new=paste0("AT.",OCID))
  
  ## Scale record ----
  if(tolower(Scale)=="relative"){
    Ao <- apply(X, 2, function(x) { max(abs(x))})
    X <-X[, .(sapply(.SD, function(x){x/max(abs(x))}))]
  }
  if(tolower(Scale)=="absolute"){
    Ao <- max(abs(X))
    X <-X[, .(sapply(.SD, function(x){x/Ao}))]
  }
  # Detrend
  X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]





  ## Resample ----
  if(Resample){
    TargetFs <- as.integer(max(2.5,kNyq)*Fmax) # 80/200 Hz
    X <- .resample(X,dt=dt,TargetFs=TargetFs,Fmax=Fmax,NW=NW,OVLP=OVLP)
    Fs <- TargetFs
    dt <- 1/Fs
  }
  # browser()
  ## Case #2. Acceleration Time Histories ----
  if(OrderTS==2){
    Wo <- X[,.(sapply(.SD, function(x) {.taperA(x,Astop=Astop.AT,Apass=Apass.AT)}))]
    AT <- X
    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(AT) <- OCID
    # browser()
    VT <- AT[, lapply(.SD, function(x){ .integrate(dx=x,dt=dt,NW=NW,OVLP=OVLP) })]

    if(nrow(VT) > nrow(Wo)){
      O <- data.table(sapply(Wo, function(x) rep(0, nrow(VT)-nrow(Wo))))
      Wo <- rbind(Wo, O)
    }
    VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
    VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(VT) <- OCID

    DT <- VT[, lapply(.SD, function(x){ .integrate(dx=x,dt=dt,NW=NW,OVLP=OVLP) })]

    if(nrow(DT) > nrow(Wo)){
      O <- data.table(sapply(Wo, function(x) rep(0, nrow(DT)-nrow(Wo))))
      Wo <- rbind(Wo, O)
    }
    DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]

    # if(0 %in% OrderEMD){
    #   DT <-DT[, .(sapply(.SD, function(x){.detrend(X=x,dt=dt,removeIMF1=removeIMF1,removeIMFn=removeIMFn)}))]
    # } else {
    #   DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]}

    DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]
    names(DT) <- OCID
  }

  # Case #1. Velocity Time Histories ----
  if(OrderTS==1){
    VT <- X
    VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]
    DT <- VT[, lapply(.SD, function(x){.integrate(dx=x,dt=dt,NW=NW,OVLP=OVLP) })]


    names(DT) <- OCID


  }

  # Case #0. Displacement Time Historyes
  if(OrderTS==0){
    DT <- X
    # Not sure if we can remove the mean in displacements.
    # DT <- DT[, .(sapply(.SD, function(x){x-mean(x)}))]

  }


  ## Derivate DT ----
  VT <- DT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
  VT <- VT[, .(sapply(.SD, function(x){x-mean(x)}))]

  ## Derivate VT ----
  AT <-  VT[, lapply(.SD, function(x){.derivate(X=x,dt=dt) })]
  AT <- AT[, .(sapply(.SD, function(x){x-mean(x)}))]

    ## Taper Zeros ----
  # Wo <- AT[,.(sapply(.SD, function(x) {.taperI(x)}))]
  Wo <- AT[,.(sapply(.SD, function(x) {.taperA(x,Astop=Astop.AT,Apass=Apass.AT)}))]
  AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]
  DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Wo[[i]]})]


  # Trim Zeros
  idx <- apply(Wo!=0,MARGIN=1,function(x){all(x)})
  AT <- AT[idx]
  VT <- VT[idx]
  DT <- DT[idx]
  # Fix trend
  AT <-AT[, .(sapply(.SD, function(x){x-mean(x)}))]
  VT <-VT[, .(sapply(.SD, function(x){x-mean(x)}))]
  DT <-DT[, .(sapply(.SD, function(x){x-mean(x)}))]


  ## Restore Scale ----
  if(tolower(Scale)=="relative"){
    AT <- AT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
    VT <- VT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
    DT <- DT[, lapply(seq_along(.SD), function(i) {.SD[[i]] * Ao[i]})]
  }
  if(tolower(Scale)=="absolute"){
    AT <- AT[, lapply(.SD, function(x) {x * Ao})]
    VT <- VT[, lapply(.SD, function(x) {x * Ao})]
    DT <- DT[, lapply(.SD, function(x) {x * Ao})]
  }
  names(AT) <- OCID
  names(VT) <- OCID
  names(DT) <- OCID

  ## Pack Time Series  ----
  NP <-  nrow(AT)
  ts <- seq(0,dt*(NP-1),dt)
  TSW <- data.table(ts=ts, AT = AT, VT = VT, DT = DT)
  ivars <- c("ts")
  mvars <- colnames(TSW[, -c("ts")])
  AUX <- data.table::melt(TSW, id.vars = ivars, measure.vars = mvars) |> na.omit()
  TSL <- AUX[,.(t=ts,s=value,ID=gsub("\\..*$", "", variable), OCID=gsub("^[^.]*\\.", "", variable))]
  ## Trim Records,  ----
  # if(TrimZeros){
  #   #Warning: Different time scales from now on
  #   TSL <- TSL[,.trimZeros(.SD),by=c("OCID","ID")]}

  ## Return ----
  Fs <- 1 / dt
  df <- Fs / NW # 0.03125#
  fs <- seq(from = 0, by = df, length.out = NW / 2)

  return(list(RTSW=RTSW,TSL=TSL,TSW=TSW,Wo=Wo,Ao=Ao,Fs = Fs, dt = dt, df=df,fs=fs,NP = NP, TargetUnits = TargetUnits, SourceUnits = Units))
}
