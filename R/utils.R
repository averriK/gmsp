#' @importFrom stats na.omit
#' @importFrom stats fft
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom seewave stdft
#' @importFrom seewave istft
#' @importFrom utils tail
#' @importFrom digest digest
#' @importFrom spectral spec.fft
#' @importFrom EMD emd
#' @importFrom hht CEEMD
#' @importFrom signal resample
#' @importFrom pracma detrend
#' @importFrom stringr str_split
#' @importFrom purrr map
#' 
#' @noRd
#'

.integrate <- function(dx,dt,Fmax=16,NW=1024,OVLP=75){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  OCID <- names(dx)
  Fpass_LP <- Fmax # 20/25 Hz
  Fstop_LP <- 1.2*Fmax # 25/30 Hz

  # Flat Zeros
  # Wo <- .taperI(dX)
  # dX <- dX*Wo

  # Build filter
  Fs <- 1/dt #
  df <- Fs / NW #
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  LP <- .buildLowPassButtterworth(f = fs, Fstop = round(1 * Fstop_LP / df) * df, Fpass = round(1 * Fpass_LP / df) * df, Astop = 0.01, Apass = 0.95)
  HI <- .buildIntegrateFilter(f = fs) ## Integrate Filter
  # Pad Zeros
  dX <- .padZeros(dx,NW=NW,OVLP=OVLP)
  # Integrate

  X <- .ffilter(dX, f = Fs, wl = NW, ovlp = OVLP, custom = HI) * NW

  return(X)
}

.resample <- function(X,dt,TargetFs,Fmax,NW=1024,OVLP=75){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  . <- NULL
  OCID <- names(X)
  Fs=1/dt

  Fpass_LP <- Fmax # 20/25 Hz
  Fstop_LP <- 1.2*Fmax # 25/30 Hz
  df <- Fs / NW #
  fs <- seq(from = 0, by = df, length.out = NW / 2)
  LP <- .buildLowPassButtterworth(f = fs, Fstop = round(1 * Fstop_LP / df) * df, Fpass = round(1 * Fpass_LP / df) * df, Astop = 0.001, Apass = 0.95)
  X <- X[, lapply(.SD, function(x) {
    ffilter(wave = x, f = Fs, wl = NW, ovlp = OVLP, custom = LP, rescale = TRUE)
  })]
  X <- X[, lapply(.SD, function(x) {
    signal::resample(x, TargetFs, Fs)
  })]
  X <-X[, .(sapply(.SD, function(x){x-mean(x)}))]
  names(X) <- OCID
  return(X)
}

.derivate <- function(X,t=NULL,dt=NULL){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  n <- length(X)
  dX <- numeric(n)
  stopifnot(n>4)
  if(!is.null(t)){
    dt <- diff(t) |> mean()
  }
  t <- seq(0,(n-1)*dt,by=dt)


  # Four-point central difference for the bulk
  for (i in 3:(n-2)) {
    dX[i] <- (-X[i+2] + 8*X[i+1] - 8*X[i-1] + X[i-2]) / (12 * dt)
  }

  # Three-point differences for the edges
  dX[1] <- (-3*X[1] + 4*X[2] - X[3]) / (2 * dt)
  dX[2] <- (-X[4] + 4*X[3] - 3*X[2]) / (2 * dt)
  dX[n-1] <- (3*X[n-1] - 4*X[n-2] + X[n-3]) / (2 * dt)
  dX[n] <- (3*X[n] - 4*X[n-1] + X[n-2]) / (2 * dt)

  return(dX)

}

.trimZeros <- function(.SD,COL="s",offset=0){
  n <- nrow(.SD)
  idx <- which(.SD[[COL]]!=0) |> first()
  START <- max(idx,2+offset)
  idx <- which(.SD[[COL]]!=0) |> last()
  END <- min(idx,n-offset)
  return(.SD[(START-offset):(END+offset)])


}

.taperA <- function(x,Astop=1e-4,Apass=1e-3){
  stopifnot(is.vector(x))
  n <- length(x)
  iH_stop <- which(abs(x) > Astop) |> first()
  iH_stop <- ifelse(is.na(iH_stop), 1, iH_stop)

  iH_pass <- which(abs(x) > Apass) |> first()
  iH_pass <- ifelse(is.na(iH_pass), 1, iH_pass)
  iH_pass <- max(iH_pass, iH_stop+1L)
  stopifnot( iH_pass > iH_stop)


  iL_stop <- which(abs(x) > Astop) |> last()
  iL_stop <- ifelse(is.na(iL_stop), n, iL_stop)

  iL_pass <- which(abs(x) > Apass) |> last()
  iL_pass <- ifelse(is.na(iL_pass), n, iL_pass)
  iL_pass <- min(iL_pass, iL_stop-1L)
  stopifnot( iL_pass < iL_stop)

  HP <- .buildHighPassButtterworth(f = seq(1, n), Fpass = iH_pass, Fstop = iH_stop, Astop = 0.01, Apass = 0.95)
  LP <- .buildLowPassButtterworth(f = seq(1, n), Fstop = iL_stop, Fpass = iL_pass, Astop = 0.01, Apass = 0.95)
  return(LP * HP)
}

.taperI <- function(x,Hstop=0.0005,Hpass=0.005,Lstop=0.999,Lpass=0.99){
  n <- length(x)
  cumIA <- cumsum((x*x)/as.numeric(x %*% x))
  iL_stop <- which(cumIA>=Lstop)|> first()#0.995
  iL_pass <- which(cumIA>=Lpass)|> first()#0.99
  LP <- .buildLowPassButtterworth(f=seq(1,n),Fstop = iL_stop, Fpass = iL_pass, Astop=0.01, Apass = 0.99)
  iH_stop <- which(cumIA>=Hstop)|> first()#0.01
  iH_pass <- which(cumIA>=Hpass)|> first()#0.005
  HP <- .buildHighPassButtterworth(f=seq(1,n),Fpass = iH_pass, Fstop = iH_stop, Astop=0.01, Apass = 0.99)
  return(HP*LP)

}

.padZeros <- function(x,nz=0,OVLP=75,NW=1024){
  on.exit(expr={rm(list = ls())}, add = TRUE)
 if(is.data.table(x)){
   NP <- nrow(x)}
  if(is.vector(x)){
     NP <- length(x)}

  NO <- OVLP*NW/100
  NB   <- ceiling((NP-NO)/(NW-NO))
  TargetNP  <- NB*(NW-NO)+NO
  NZ <- max(nz,TargetNP-NP)
  if(is.data.table(x) & length(NZ)>0){
    O <- data.table(sapply(x, function(x) rep(0, NZ)))
    x <- rbind(x, O)
  }

  if(is.vector(x) & length(NZ)>0){
    O <- rep(0, NZ)
    x <- c(x, O)
  }
  return(x)


}

.ffilter <- function (x, f, custom, ovlp=75,wl=1024) {
  on.exit(expr={rm(list = ls())}, add = TRUE)
  NP <- length(x)
  # browser()
  STEP <- seq(1, NP + 1 - wl, wl - (ovlp * wl/100))
  z <- (seewave::stdft(
    wave = matrix(x,ncol=1), f = f, wl = wl, zp = 0, step = STEP,
    wn = "hanning", fftw = FALSE, complex = TRUE))
  z <- z*custom
  X <- (seewave::istft(z, wl = wl, ovlp = ovlp, wn = "hanning", output = "matrix",f = f)) |> as.vector()
  return(X)
}

.buildLowPassButtterworth <- function(f=NULL,Fstop=NULL,Fpass=NULL,Astop=NULL,Apass=NULL) {
  on.exit(expr={rm(list = ls())}, add = TRUE)
  OK <- !is.null(f)&& is.vector(f) && !is.null(Fstop) && !is.null(Fpass) && Fpass<Fstop && !is.null(Astop) && !is.null(Apass) && Astop<Apass
  stopifnot(OK)
  TOL <- Astop
  N <- .NminLP(fS=Fstop,fP=Fpass,aS=Astop,aP=Apass)
  LP <- 1/sqrt(1-(Astop^2-1)*((f/Fstop)^(2*N))/(Astop^2))
  LP[LP<TOL] <- 0
  # NminHP[fs, fp, As, Ap] == NminLP[fs, fp, 1 - As, 1 - Ap]
  return(LP)
}

.buildHighPassButtterworth <- function(f=NULL,Fstop=NULL,Fpass=NULL,Astop=NULL,Apass=NULL) {
  on.exit(expr={rm(list = ls())}, add = TRUE)

  OK <- !is.null(f)&& is.vector(f) && !is.null(Fstop) && !is.null(Fpass) && Fpass>Fstop && !is.null(Astop) && !is.null(Apass) && Astop<Apass
  stopifnot(OK)
  TOL <- Astop
  n <- .NminLP(fS=Fstop,fP=Fpass,aS=1-Astop,aP=1-Apass)
  HP <- 1 - 1/sqrt(1 + (-1 + (-1 + Astop)^(-2))*(f/Fstop)^(2*n))
  HP[HP<TOL] <- 0
  # NminHP[fs, fp, As, Ap] == NminLP[fs, fp, 1 - As, 1 - Ap]
  return(HP)
}

.NminLP <- function(fS=0,fP=0,aS=NULL,aP=NULL) {
  on.exit(expr={rm(list = ls())}, add = TRUE)
  OK <-     fS>0 && fP>0 && !is.null(aS) && !is.null(aP)
  stopifnot(OK)
  n <- log((aS*sqrt((-1 + aP^2)/(-1 + aS^2)))/aP)/log(fP/fS)

  return(n)
}

.buildIntegrateFilter <- function(f){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  ws  <- 2*pi*f
  HI <- c(0,1/(1i*ws[2:length(ws)]))
  # HI <- purrr::prepend(values=0,1/(1i*ws[2:length(fs)]))
  return(HI)
}

.buildDerivateFilter <- function(f){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  ws  <- 2*pi*f
  HD <- c(0,1i*ws[2:length(ws)])
  # HD <- purrr::prepend(values=0,(1i*ws[2:length(fs)]))
  return(HD)
}

.getSF <- function(SourceUnits,TargetUnits,g_mms2=9806.650){
  switch (SourceUnits,
          "m" = switch(TargetUnits,"mm"= 1000,"cm"=100,"m"=1,NULL),
          "cm" = switch(TargetUnits,"mm"= 10,"cm"=1,"m"=1/100,NULL),
          "mm" = switch(TargetUnits,"mm"= 1,"cm"=1/10,"m"=1/1000,NULL ),
          "g" = switch (TargetUnits,"m" = g_mms2/1000,"cm" = g_mms2/10,"mm" = g_mms2/1,NULL),
          "gal" = switch(TargetUnits,"mm"= 10,"cm"=1,"m"=1/100,NULL),
          NULL
  )
}

.getG <- function(TargetUnits,g_mms2=9806.650){
  switch (TargetUnits,
          "m" = g_mms2/1000,
          "cm" = g_mms2/10,
          "mm" = g_mms2/1,
          NULL
  )
}







.getFFT <- function(.SD){
  FFT <- spectral::spec.fft(y=.SD$Y,x=.SD$X,center = TRUE)
  data.table(f=FFT$fx,PSD=(FFT$PSD))
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


.isRational <- function(x) {
  # Check if the number is finite
  if (!is.finite(x)) {
    return(FALSE)
  }
  
  # Convert the number to a character string
  XS <- as.character(x)
  
  # Check for infinite decimal expansion
  if (grepl("\\.\\d*", XS)) {
    # Look for a repeating decimal (like 0.33333...)
    if (any(grepl("(\\d)\\1{2,}", XS))) {
      return(FALSE)
    }
    
    # If it does not have a repeating decimal, it is rational
    return(TRUE)
  }
  
  # If no decimal part is present, it is rational (like an integer)
  return(TRUE)
}

# Examples
# is_rational(0.5)    # TRUE
# is_rational(0.3333) # TRUE
# is_rational(1/3)    # FALSE (not expressible in finite decimal form)
# is_rational(5)      # TRUE (integer, hence rational)
# is_rational(NaN)    # FALSE
# is_rational(Infinity) # FALSE



