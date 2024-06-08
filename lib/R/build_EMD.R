
#' Title
#'
#' @param s numeric vector
#' @param t numeric vector
#' @param dt numeric
#' @param method string
#' @param boundary string
#' @param max.imf integer
#' @param noise.type string
#' @param noise.amp numeric
#' @param trials integer
#' @param stop.rule string
#' @param plot boolean
#' @param verbose boolean
#'
#' @importFrom data.table :=
#' @importFrom stats na.omit
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom data.table melt
#' @importFrom hht EEMDCompile
#' @importFrom hht CEEMD
#' @importFrom hht EEMD
#' @importFrom EMD emd
#'
#' @return list
#' @export 
#'
#' @examples
#'
build_EMD <- function(s,t=NULL,dt=NULL,method="emd",boundary="wave", max.imf=15,noise.type="gaussian",noise.amp=0.5e-7,trials=10,stop.rule="type5",plot=TRUE,verbose=FALSE){
  on.exit(expr = {rm(list = ls())}, add = TRUE)
stopifnot(!is.null(s) && tolower(method) %in% c("emd","eemd","ceemd") && tolower(stop.rule) %in% c("type1","type2","type3","type4","type5") && tolower(boundary) %in% c("none","wave","symmetric","periodic","evenodd") && noise.type %in% c("uniform","gaussian") )

  . <- NULL
  n <- length(s)
  stopifnot(n>4)
  if(!is.null(t)){
    dt <- diff(t) |> mean()
  }
  t <- seq(0,(n-1)*dt,by=dt)
  # Trim Zeros
  DT <- data.table(t=t,s=s)
  DT <- .trimZeros(DT)
  if(tolower(method)=="emd"){
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf,stoprule=stop.rule)
    IMF <- AUX$imf  |> as.data.table()
    names(IMF) <- paste0("IMF", seq_len(ncol(IMF)))
    RES <- AUX$residue|> as.vector() |> unname()

  }

  if(tolower(method)=="eemd"){
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf,stoprule=stop.rule)
    nimf <- AUX$nimf
    DIR <- tempdir(check = TRUE)
   hht::EEMD(sig=DT$s, tt=DT$t,nimf=nimf,max.imf=max.imf,boundary=boundary,noise.amp=noise.amp, noise.type=noise.type,trials=trials,stop.rule=stop.rule,trials.dir = DIR,verbose=verbose)
   AUX <- EEMDCompile(trials.dir = DIR, trials=trials, nimf=nimf) |> suppressWarnings()
   unlink(DIR,force = TRUE,recursive = TRUE)
   IMF <- AUX$averaged.imfs |> as.data.table()
   names(IMF) <- paste0("IMF", seq_len(ncol(IMF)))
    RES <- AUX$averaged.residue|> as.vector() |> unname()

  }

  if(tolower(method)=="ceemd"){
    AUX <- hht::CEEMD(sig=DT$s, tt=DT$t,noise.amp=noise.amp, noise.type=noise.type,trials=trials,stop.rule=stop.rule,verbose=verbose)
    IMF <- AUX$imf |> as.data.table()
    names(IMF) <- paste0("IMF", seq_len(ncol(IMF)))
    RES <- AUX$residue|> as.vector() |> unname()

  }


  Tm <- IMF[,lapply(.SD, function(x) {.getTm(x,Fs=1/dt)})]
  wm <- 2*pi/Tm
  fm <- 1/Tm
  PGA <- IMF[,lapply(.SD, function(x) {max(abs(x))})]
  nimf <- ncol(IMF)
  DATA <- NULL
  if(plot==TRUE){

    M <- data.table(signal=DT$s,IMF,residue=RES)
    offset <- 1.25*ceiling(max(M)-min(M))
    for(i in 1:ncol(M)){
      j <- ncol(M)-i+1
      M[[j]] <- M[[j]]+offset*i
    }

    AUX <- data.table(t=DT$t,M)
    ivars <- c("t")
    mvars <- colnames(AUX[, -c("t")])
    DATA <- melt(AUX, id.vars = ivars, measure.vars = mvars) |> na.omit()
    DATA <- DATA[,.(X=t,Y=value,ID=variable)]
    rm(M)
  }

  # Restore zeros
  nimf <- ncol(IMF)
  i <- as.integer(DT$t[1]/dt)
  j <- as.integer(last(DT$t)/dt)

  OB <- data.table(matrix(0, nrow = n-j, ncol = ncol(IMF)))
  OA <- data.table(matrix(0, nrow = i-1, ncol = ncol(IMF)))
  IMF <- rbindlist(list(OA,IMF,OB),use.names = FALSE)


  return(list(s=s,t=t,fm=fm,Tm=Tm,pga=PGA,imf=IMF,nimf=nimf,residue=RES,plot.data=DATA,plot.offset=offset))
}

