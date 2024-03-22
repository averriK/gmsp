
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
#'
#' @return list
#' @export buildIMF
#'
#' @examples
#'
buildIMF <- function(s,t=NULL,dt=NULL,method="emd",boundary="wave", max.imf=20,noise.type="gaussian",noise.amp=0.5e-7,trials=10,stop.rule="type5",plot=TRUE){
  on.exit(expr = {rm(list = ls())}, add = TRUE)
stopifnot(!is.null(s) && tolower(method) %in% c("emd","semd","eemd","ceemd") && tolower(stop.rule) %in% c("type1","type2","type3","type4","type5") && tolower(boundary) %in% c("wave","spline","mirror","extrapolate") && noise.type %in% c("uniform","gaussian") )

  . <- NULL
  if(is.null(dt) & !is.null(t)){
    dt <- t[2]-t[1]
  }

  if(is.null(t) & !is.null(dt)){
    n <- length(s)
    t <- seq(0,(n-1)*dt,by=dt)
  }

  stopifnot(dt == t[2]-t[1] & length(s)==length(t) )
  # Trim Zeros

  DT <- data.table(t=t,s=s)
  n <- nrow(DT)
  DT <- .trimZeros(DT)

  if(tolower(method)=="emd"){
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf,stoprule=stop.rule)
    M <- AUX$imf  |> as.data.table()
    NC <- ncol(M)
    RES <- M[[NC]] |> unname()
    IMF <- M[,-NC,with = FALSE]
  }

  if(tolower(method)=="eemd"){
    AUX <- EMD::emd(xt=DT$s, tt=DT$t, boundary=boundary, max.imf=max.imf,stoprule=stop.rule)
    nimf <- AUX$nimf
    DIR <- tempdir()
   hht::EEMD(sig=DT$s, tt=DT$t,nimf=nimf,max.imf=max.imf,boundary=boundary,noise.amp=noise.amp, noise.type=noise.type,trials=trials,stop.rule=stop.rule,trials.dir = DIR)
   AUX <- EEMDCompile(trials.dir = DIR, trials=trials, nimf=nimf) |> suppressWarnings()
   unlink(DIR,force = TRUE,recursive = TRUE)
   IMF <- AUX$averaged.imfs |> as.data.table()
    RES <- AUX$averaged.residue |> unname()
    M <- data.table(IMF,RES)
  }

  if(tolower(method)=="ceemd"){
    AUX <- hht::CEEMD(sig=DT$s, tt=DT$t,noise.amp=noise.amp, noise.type=noise.type,trials=trials,stop.rule=stop.rule)
    IMF <- AUX$imf |> as.data.table()
    RES <- AUX$residue |> unname()
    M <- data.table(IMF,RES)
  }
  names(IMF) <- paste0("IMF", seq_len(ncol(IMF)))

  # Tm <- sapply(IMF, function(x) {.getTm(x,Fs=1/dt)})
  Tm <- IMF[,lapply(.SD, function(x) {.getTm(x,Fs=1/dt)})]
  wm <- 2*pi/Tm
  fm <- 1/Tm
  PGA <- IMF[,lapply(.SD, function(x) {max(abs(x))})]
  nimf <- ncol(IMF)
  DATA <- NULL
  if(plot==TRUE){
    M <- IMF
    offset <- 1.25*ceiling(max(M)-min(M))

    for(i in 1:nimf){
      j <- nimf-i+1
      M[[j]] <- M[[j]]+offset*i
    }

    AUX <- data.table(t=DT$t,"Residue"=RES,"Signal"=DT$s+offset*(nimf+2),M)
    ivars <- c("t")
    mvars <- colnames(AUX[, -c("t")])
    DATA <- melt(AUX, id.vars = ivars, measure.vars = mvars) |> na.omit()
    DATA <- DATA[,.(X=t,Y=value,ID=variable)]
  }
  # Restore zeros
  nimf <- ncol(IMF)
  i <- as.integer(DT$t[1]/dt)
  j <- as.integer(last(DT$t)/dt)

  RES=c(rep(0,times=i-1),RES)
  RES <- c(RES,rep(0,times=n-length(RES)))

  Oi <- data.table(matrix(0, nrow = n-j, ncol = ncol(IMF)))
  Oj <- data.table(matrix(0, nrow = i-1, ncol = ncol(IMF)))
  IMF <- rbindlist(list(Oj,IMF,Oi),use.names = FALSE)
  names(IMF) <- paste0("IMF", seq_len(ncol(IMF)))

  return(list(s=s,t=t,fm=fm,Tm=Tm,pga=PGA,imf=IMF,nimf=nimf,residue=RES,plot.data=DATA,plot.offset=offset))
}

