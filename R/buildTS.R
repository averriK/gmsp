
#' Title
#'
#' @param a data.table Time Series
#' @param dt numeric Time Step
#' @param UN character Units
#' @param Amin numeric Minimum Amplitude
#' @param Amax numeric Maximum Amplitude
#' @param TFT data.table Transfer Functions Table
#' @param DownFs integer Downsample Frequency
#' @param UpFs integer Upsample Frequency
#' @param FlatZerosAT boolean Flat Zeros Acceleration Time Series
#' @param FlatZerosVT boolean Flat Zeros Velicity Time Series
#' @param FlatZerosDT boolean Flat Zeros Displacements Time Series
#' @param DerivateDT boolean Derivate Displacements Time Series
#' @param DerivateVT boolean Derivate Velocity Time Series
#' @param DetrendAT boolean Detrend Acceleration Time Series
#' @param DetrendVT boolean Detrend Velocity Time Series
#' @param DetrendDT boolean Detrend Displacement Time Series
#' @param Fpass_LP numeric Low Pass Frequency
#' @param Fstop_LP numeric Low Stop Frequency
#' @param Fpass_HP numeric High Pass Frequency
#' @param Fstop_HP numeric High Stop Frequency
#' @param RestoreScale boolean
#' @param TargetUnits character Units
#' @param NW integer Windows Length
#' @param OVLP integer Overlap
#'
#' @return data.table
#' @export buildTS
#'
#' @examples
#'
#'
#'
#' @import data.table
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom seewave ffilter
#' @importFrom signal resample
#' @importFrom pracma detrend
#' @importFrom stringr str_split
#' @importFrom purrr map
#'
#'
buildTS <- function(a,dt,UN=NULL,Amin = 0,Amax=Inf,TFT=NULL,DownFs=0,UpFs=0,FlatZerosAT=FALSE,FlatZerosVT=FALSE,FlatZerosDT=FALSE,DerivateDT=FALSE,DerivateVT=FALSE,DetrendAT=TRUE,DetrendVT=FALSE,DetrendDT=FALSE,Fpass_LP=0,Fstop_LP=0,Fpass_HP=0,Fstop_HP=0,RestoreScale=FALSE,TargetUnits="mm",NW=2048,OVLP=75){
  on.exit(expr={rm(list = ls())}, add = TRUE)


  OK <- is.data.table(a) #&& !is.null(dt) && !is.null(UN)
  stopifnot(OK)

  ATo <- copy(a)


  ## Check Record ----------------------------------------------------------------------
  if(
    nrow(ATo)==0||
    dt==0||
    any(is.na(ATo))||
    max(abs(ATo))==0 #||
    #any(sapply(ATo, function(x){max(abs(x))})==0)
  ){
    # Null Record
    return(NULL)
  }

  OCID <- names(ATo)
  if(any(grepl(OCID,pattern = "^T$"))){
    ATo[,c("T"):=NULL]
    OCID <- names(ATo)
  }
  if(length(unique(OCID))<3){
    # Invalid Header
    return(NULL)
  }
  if(grepl(UN,pattern = "[///+]")){
    UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
  }
  if(tolower(UN) %in% c("mm","cm","m","gal","g")){
    ## Scale Units ----------------------------------------------------------------------
    if(UN!=TargetUnits) {
      SFU <- .getSF(SourceUnits =tolower(UN),TargetUnits = TargetUnits)
      # AT <- map(AT,function(x){x*SFU})
      ATo[,(colnames(ATo)):=lapply(.SD,function(x){x*SFU})]
    } else {
      SFU <- 1
    }
  }

  Fs <- 1/dt

  ## Set Scale Reference ----------------------------------------------------------------------
  DUMMY <- NULL
  PGAo <- apply(ATo,2,function(x){max(abs(x))})


  # Keep Only Records with PGA > Amin & PGA < Amax ------
  if(max(PGAo)<Amin | max(PGAo)>Amax) return(NULL)

  ## Check Length ----------------------------------------------------------------------
  NP <- nrow(ATo)
  if(NP<NW){return(NULL)}

  ## Detrend ----------------------------------------------------------------------
  if(DetrendAT){
    ATo[,(colnames(ATo)):=lapply(.SD,function(x){pracma::detrend(x,tt="linear")})]
  }


  ## Flat Zeros & Taper  ---------------------------------------------------------------------
  if(FlatZerosAT==TRUE){
    ATo <- ATo[,lapply(X=.SD,FUN= function(x){
      n <- length(x)
      iH_stop <- which(abs(x)>0.1) |> first() # 0.1 mm/s2
      iL_stop <- which(abs(x)>0.1) |> last()
      iH_pass <- which(abs(x)>1) |> first() #1 mm/s2
      iL_pass <- which(abs(x)>1) |> last()

      if(length(iH_pass)==1 && length(iH_stop)==1 && iH_pass>iH_stop){
        HP <- .buildHighPassButtterworth(f=seq(1,n),Fpass = iH_pass, Fstop = iH_stop, Astop=0.001, Apass = 0.999)
      } else {
        HP <- rep(1,n)}
      if(length(iL_pass)==1 && length(iL_stop)==1 && iL_pass<iL_stop){
        LP <- .buildLowPassButtterworth(f=seq(1,n),Fstop = iL_stop, Fpass = iL_pass, Astop=0.001, Apass = 0.999)
      } else {
        LP <- rep(1,n)}
      return(x*LP*HP)
    })]
  }



  ## Up-sampling ----------------------------------------------------------------------

  if(UpFs!=0 & UpFs>Fs){
    ATo <- ATo[,lapply(.SD,function(x){
      x <- signal::resample(x,UpFs,Fs)
      # x <- pracma::detrend(x,tt="linear")
      return(x)
    })]

    ## Detrend ----------------------------------------------------------------------
    if(DetrendAT){
      ATo[,(colnames(ATo)):=lapply(.SD,function(x){pracma::detrend(x,tt="linear")})]
    }

    # Update Record
    Fs <- UpFs
    dt <- 1/Fs
    NP <-  nrow(ATo)
    names(ATo) <- OCID
  }

  ## Padding Zeros ----------------------------------------------------------------------
  NP <- nrow(ATo)
  NZ <- .getNZ(NP)
  if(NZ > 0) {
    ZEROS <- data.table()[,(colnames(ATo)):=list(rep(0,NZ))]
  } else {
    ZEROS <- data.table()[,(colnames(ATo)):=list(rep(0,NW))]
  }
  ATo <- rbindlist(list(ZEROS,ATo))



  ## Antialias & Downsampling ----------------------------------------------------------------------

  if(DownFs!=0 & Fs!=DownFs){
    df <- Fs/NW #0.03125#
    fs <- seq(from=0,by=df,length.out=NW/2)
    FsNYQ <- Fs/2 #16
    COLS <- colnames(ATo)

    if(Fstop_LP==0 & Fpass_LP==0){
      LP <- rep(1,times=NW/2)}
    else {
      LP  <- .buildLowPassButtterworth(f=fs,Fstop = round(1*Fstop_LP/df)*df, Fpass=round(1*Fpass_LP/df)*df,Astop = 0.001,Apass = 0.95)}


    ATo <- ATo[,lapply(.SD,function(x){
      x <- ffilter(wave=x, f = Fs, wl = NW, ovlp = OVLP,custom = LP,rescale = TRUE)
      x <- signal::resample(x,DownFs,Fs)
      return(x)
    })]

    ## Detrend ----------------------------------------------------------------------
    if(DetrendAT){
      ATo <- ATo[,lapply(.SD,function(x){
        x <- pracma::detrend(x,tt="linear")
        return(x)
      })]
    }
    names(ATo) <- COLS
    # Update
    Fs <- DownFs
    dt <- 1/Fs

  }

  ## Padding zeros ----------------------------------------------------------------------
  NP <- nrow(ATo)
  if(NP<NW){return(NULL)}
  NZ <- .getNZ(NP)
  if(NZ > 0) {
    ZEROS <- data.table()[,(colnames(ATo)):=list(rep(0,NZ))]
  } else {
    ZEROS <- data.table()[,(colnames(ATo)):=list(rep(0,NW))]
  }
  ATo <- rbindlist(list(ZEROS,ATo))

  ## Build HT Filters ----------------------------------------------------------------------
  Fs <- 1/dt
  df <- Fs/NW #0.03125#
  fs <-  seq(from=0,by=df,length.out=NW/2)

  # SiteSN <- colnames(TFT)


  ## Convolute ----------------------------------------------------------------------
  # COLS <- colnames(TFT)#c("I",names(TFT)|> grep(pattern = "HS2O",value = TRUE))
  if(!is.null(TFT)){
    # TFT <- data.table(I=rep(1,times=nrow(ATo)))
    AT <- map(TFT,function(x){
      map(ATo,function(y){
        y <- .ffilter(y, f = Fs, wl = NW, ovlp = OVLP, custom = x)*NW
        # y <- pracma::detrend(y,tt="linear")
        return(y)
      }) |> as.data.table()
    })
    names(AT) <- colnames(TFT)
    ## GROUP TIMESERIES IN DT
    AT <- as.data.table(AT)
  } else {
    AT <- as.data.table(ATo)
  }

  ## Restore Scale -------------------------------------------------------------------

  if(RestoreScale){
    ICOLS <- (colnames(AT) |> grep(pattern = "^I.",value = TRUE))
    if(length(ICOLS)>0){
      PGA <- apply(AT[,..ICOLS],2,function(x){max(abs(x))})
    } else {
      PGA <- apply(AT,2,function(x){max(abs(x))})
    }
    SF <- (PGAo/PGA)
    for(n in seq_len(length.out = length(SF))){
      COLS <- (colnames(AT) |> grep(pattern = OCID[n],value = TRUE))
      AT[,(COLS):=lapply(.SD,function(x){x*SF[n]}),.SDcols=COLS]
    }
  }




  ## Padding zeros ----------------------------------------------------------------------
  NP <- nrow(AT)
  NZ <- .getNZ(NP)
  if(NZ > 0) {
    ZEROS <- data.table()[,(colnames(AT)):=list(rep(0,NZ))]
  } else {
    ZEROS <- data.table()[,(colnames(AT)):=list(rep(0,NW))]
  }
  AT <- rbindlist(list(ZEROS,AT))
  ## Build Filters ----------------------------------------------------------------------
  Fs <- 1/dt
  df <- Fs/NW #0.03125#
  fs <-  seq(from=0,by=df,length.out=NW/2)


  FsNYQ <- Fs/2 #30
  if(Fstop_LP==0 & Fpass_LP==0){
    LP <- rep(1,times=NW/2)} else{
      LP <-  .buildLowPassButtterworth(f=fs, Fstop = round(1*Fstop_LP/df)*df, Fpass=round(1*Fpass_LP/df)*df, Astop=0.001, Apass=0.90)
    }

  #
  if(Fpass_HP==0 & Fstop_HP==0){
    HP <- rep(1,times=NW/2)} else {
      HP <-  .buildHighPassButtterworth(f=fs, Fstop = round(1*Fstop_HP/df)*df, Fpass=round(1*Fpass_HP/df)*df,Astop=0.001, Apass=0.99)
    }
  ## Highpass & Lowpass filters`------------------------------------------------------------------
  COLS <- colnames(AT)
  AT <- AT[,lapply(.SD,function(x){
    x <- ffilter(wave=x, f = Fs, wl = NW, ovlp = OVLP,custom = HP*LP,rescale = TRUE)
    # x <- pracma::detrend(x,tt="linear")
    return(x)
  })]
  ## Detrend ----------------------------------------------------------------------
  if(DetrendAT){
    AT <- AT[,lapply(.SD,function(x){
      x <- pracma::detrend(x,tt="linear")
      return(x)
    })]
  }


  names(AT) <- COLS

  ## Restore Scale -------------------------------------------------------------------

  if(RestoreScale){
    ICOLS <- (colnames(AT) |> grep(pattern = "^I.",value = TRUE))
    # PGA <- apply(AT[,..ICOLS],2,function(x){max(abs(x))})

    if(length(ICOLS)>0){
      PGA <- apply(AT[,..ICOLS],2,function(x){max(abs(x))})
    } else {
      PGA <- apply(AT,2,function(x){max(abs(x))})
    }

    SF <- (PGAo/PGA)
    for(n in seq_len(length.out = length(SF))){
      COLS <- (colnames(AT) |> grep(pattern = OCID[n],value = TRUE))
      AT[,(COLS):=lapply(.SD,function(x){x*SF[n]}),.SDcols=COLS]
    }}




  ## Build VT Filters ----------------------------------------------------------------------
  Fs <- 1/dt
  df <- Fs/NW #0.03125#
  fs <-  seq(from=0,by=df,length.out=NW/2)

  FsNYQ <- Fs/2 #16
  if(Fstop_LP==0 & Fpass_LP==0){
    LP <- rep(1,times=NW/2)} else {
      LP <-  .buildLowPassButtterworth(f=fs, Fstop = round(Fstop_LP/df)*df, Fpass=round(Fpass_LP/df)*df, Astop=0.001, Apass=0.90)
    }

  HI <-  .buildIntegrateFilter(f=fs) ## Integrate Filter

  ## Padding zeros AT----------------------------------------------------------------------
  NP <- nrow(AT)
  NZ <- .getNZ(NP)
  if(NZ > 0) {
    ZEROS <- data.table()[,(colnames(AT)):=list(rep(0,NZ))]
  } else {
    ZEROS <- data.table()[,(colnames(AT)):=list(rep(0,NW))]
  }
  AT <- rbindlist(list(ZEROS,AT))
  # Integrate Acceleration ----------------------------------------------------------------------
  COLS <- colnames(AT)
  VT <- AT[,lapply(.SD,function(x){
    x <- .ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)*NW
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
    # x <- pracma::detrend(x,tt="linear")
  })]
  names(VT) <- COLS

  ## Build Displacement Filters ----------------------------------------------------------------------
  Fs <- 1/dt
  df <- Fs/NW #0.03125#
  fs <-  seq(from=0,by=df,length.out=NW/2)
  FsNYQ <- Fs/2 #16 round(3/8*Fpass_LP/df)*df
  if(Fstop_LP==0 & Fpass_LP==0){
    LP <- rep(1,times=NW/2)} else {
      LP <-  .buildLowPassButtterworth(f=fs, Fstop = round(Fstop_LP/df)*df, Fpass=round(Fpass_LP/df)*df, Astop=0.001, Apass=0.90)
    }
  HI <-  .buildIntegrateFilter(f=fs) ## Integrate Filter

  ## Padding zeros VT ----------------------------------------------------------------------
  NP <- nrow(VT)
  NZ <- .getNZ(NP)
  if(NZ > 0) {
    ZEROS <- data.table()[,(colnames(VT)):=list(rep(0,NZ))]
  } else {
    ZEROS <- data.table()[,(colnames(VT)):=list(rep(0,NW))]
  }
  VT <- rbindlist(list(ZEROS,VT))


  ## Integrate Velocity ----------------------------------------------------------------------
  COLS <- colnames(AT)
  DT <- VT[,lapply(.SD,function(x){
    x <- NW*.ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HI)
    x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
    # x <- pracma::detrend(x,tt="linear")
  })]
  names(DT) <- COLS

  ## Build Derivative Filters ----------------------------------------------------------------------
  FsNYQ <- Fs/2 #30

  if(Fstop_LP==0 & Fpass_LP==0){
    LP <- rep(1,times=NW/2)} else {
      LP <-  .buildLowPassButtterworth(f=fs, Fstop = round(1*Fstop_LP/df)*df, Fpass=round(1*Fpass_LP/df)*df, Astop=0.001, Apass=0.90)
    }

  HD <-  .buildDerivateFilter(f=fs)## Derivate Filter

  if(Fpass_HP==0 & Fstop_HP==0){
    HP <- rep(1,times=NW/2)} else{
      HP <-  .buildHighPassButtterworth(f=fs, Fstop = round(1*Fstop_HP/df)*df, Fpass=round(1*Fpass_HP/df)*df ,Astop=0.001, Apass=0.99)
    }

  # HP <-  .buildHighPassButtterworth(f=fs, Fstop = Fstop_HP,Fpass=Fpass_HP ,Astop=0.001, Apass=0.99)



  # Derivate DT   ----------------------------------------------------------------------
  if(DerivateDT){
    if(Fstop_LP==0 & Fpass_LP==0){LP <- rep(1,times=NW/2)} else {
      LP <-  .buildLowPassButtterworth(f=fs, Fstop = round(Fstop_LP/df)*df, Fpass=round(Fpass_LP/df)*df, Astop=0.001, Apass=0.90)
    }

    COLS <- colnames(AT)
    VT <- DT[,lapply(.SD,function(x){
      x <- NW*.ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
      # x <- pracma::detrend(x,tt="linear")
      return(x)
    })]
    names(VT) <- COLS
  }


  # Derivate VT   ----------------------------------------------------------------------
  if(DerivateVT){
    COLS <- colnames(AT)
    AT <- VT[,lapply(.SD,function(x){
      x <- NW*.ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HD)
      x <- seewave::ffilter(x, f = Fs, wl = NW, ovlp = OVLP, custom = HP*LP, rescale = TRUE)
      # x <- pracma::detrend(x,tt="linear")
      return(x)
    })]
    names(AT) <- COLS
  }


  # Restore Scale -------------------------------------------------------------------
  if(RestoreScale){
    ICOLS <- (colnames(AT) |> grep(pattern = "^I.",value = TRUE))
    # PGA <- apply(AT[,..ICOLS],2,function(x){max(abs(x))})
    if(length(ICOLS)>0){
      PGA <- apply(AT[,..ICOLS],2,function(x){max(abs(x))})
    } else {
      PGA <- apply(AT,2,function(x){max(abs(x))})
    }
    SF <- (PGAo/PGA)
    for(n in seq_len(length.out = length(SF))){
      COLS <- (colnames(AT) |> grep(pattern = OCID[n],value = TRUE))
      AT[,(COLS):=lapply(.SD,function(x){x*SF[n]}),.SDcols=COLS]
      VT[,(COLS):=lapply(.SD,function(x){x*SF[n]}),.SDcols=COLS]
      DT[,(COLS):=lapply(.SD,function(x){x*SF[n]}),.SDcols=COLS]
    }
  }


  ## Pack & Taper  -------------------------------------------------

  ICOLS <- colnames(AT) |> grep(pattern = "^I.",value = TRUE)

  if(length(ICOLS)>0){
    ATo <- AT[,..ICOLS] |> as.data.table() |> setnames(new = OCID)
    VTo <- VT[,..ICOLS] |> as.data.table() |> setnames(new = OCID)
    DTo <- DT[,..ICOLS] |> as.data.table() |> setnames(new = OCID)
  } else {
    ATo <- AT |> as.data.table() |> setnames(new = OCID)
    VTo <- VT |> as.data.table() |> setnames(new = OCID)
    DTo <- DT |> as.data.table() |> setnames(new = OCID)
  }

  PGA <- apply(ATo,2,function(x){max(abs(x))})

  NMX <- min(nrow(ATo),nrow(VTo),nrow(DTo))
  ATo <- ATo[-((NMX):.N)]
  VTo <- VTo[-((NMX):.N)]
  DTo <- DTo[-((NMX):.N)]
  Wo <- .taper(ATo[[which.max(PGA)]])
  ATo[,(OCID):=lapply(.SD,function(x){Wo*x}),.SDcols=OCID]
  VTo[,(OCID):=lapply(.SD,function(x){Wo*x}),.SDcols=OCID]
  DTo[,(OCID):=lapply(.SD,function(x){Wo*x}),.SDcols=OCID]
  TS <- list(data.table(W=Wo,AT=ATo,VT=VTo,DT=DTo))
  names(TS) <- "I"

  if(!is.null(TFT)){
    XCOLS <- colnames(TFT)[colnames(TFT)!="I"]#SiteSN[SiteSN!="I"]
    for(j in seq_along(along.with = XCOLS)){
      COLS <- grep(colnames(AT), pattern = paste0(XCOLS[j],"\\."), value = TRUE)
      ATo <- AT[,..COLS] |> as.data.table() |> setnames(new = OCID)
      PGA <- apply(ATo,2,function(x){max(abs(x))})
      Wo <- .taper(ATo[[which.max(PGA)]])
      VTo <- VT[,..COLS] |> as.data.table() |> setnames(new = OCID)
      DTo <- DT[,..COLS] |> as.data.table() |> setnames(new = OCID)

      ATo[,(OCID):=lapply(.SD,function(x){Wo*x}),.SDcols=OCID]
      VTo[,(OCID):=lapply(.SD,function(x){Wo*x}),.SDcols=OCID]
      DTo[,(OCID):=lapply(.SD,function(x){Wo*x}),.SDcols=OCID]
      TSo <- list(data.table(W=Wo,AT=ATo,VT=VTo,DT=DTo))
      names(TSo) <- XCOLS[j]
      TS <- append(TS,TSo)
    }
    rm(AT,ATo,VT,VTo,DT,DTo)
  }




  ## Trim Zeros ---------------------------------------------------------------------
  TS <- lapply(TS,unique)
  ## Add row of Zeros
  TS$I <- rbindlist(list(TS$I,0*TS$I[1]))

  ## Return ---------------------------------------------------------------------
  return(list(TS=TS,Fs=Fs,dt=dt,NP=NP,UN="mm"))

}
