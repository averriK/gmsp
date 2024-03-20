devtools::load_all()
library(xplot)
library(data.table)
if(!exists("SET")){
  RecordsFolder <- file.path("/Users/averri/Database/gmdb/source/tables")
  SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
}

RSN.TARGET <- 1540  #577 300 1500 1540
ID_TARGET <- "DT"
OCID_TARGET <- "H1"
COMPLETE <- FALSE

# -----
RAW <- SET[[RSN.TARGET]] #577 300 1500
RECORD <- buildTS(
  RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  Fmax=20,
  TrimZeros=FALSE,
  Detrend=TRUE,
  PadZeros=TRUE,
  TargetUnits="mm",
  NW=1024,
  OVLP=75,
  Astop=0.5e-4,
  Apass=1e-3)
TSL <- RECORD$TSL
dt <- RECORD$dt


DATA.TS <- TSL[ID=="AT" & OCID=="H1",.(X=t,Y=s,ID=paste0(ID,".",OCID))]
plot.ggplot2(DATA.TS, plot.type = "line",line.size=0.5)

# -----
s <- DATA$Y
t <- DATA$X
NIT <- 10
NIMF <- 8
n <- 6# n<0
NAMP <- 0.5*10^(-n)
NTYPE <- "gaussian" #c("uniform","gaussian")

if(COMPLETE){
  IMF <- hht::CEEMD(sig=s, tt=t, noise.amp=NAMP, trials=NIT, verbose = TRUE, noise.type=NTYPE)
  M <- IMF$imf
} else {
  TMP <- tempdir(check = TRUE)
  hht::EEMD(sig=s, tt=t, noise.amp=NAMP, trials=NIT, nimf=NIMF, trials.dir = TMP )
  #Compile the results
  IMF <- hht::EEMDCompile(trials=NIT, nimf=NIMF, trials.dir = TMP)
  unlink(TMP)
  M <- IMF$averaged.imfs
}

OFFSET <- 1.25*ceiling(max(M)-min(M))
for(i in 1:ncol(M)){
  j <- ncol(M)-i+1
  M[,j] <- M[,j]+OFFSET*i
}
AUX <- as.data.table(M)
AUX[,`:=`(t=IMF$tt,"Residue"=IMF$residue,"Signal"=IMF$original.signal+OFFSET*(ncol(M)+2))]
setnames(AUX,old=colnames(AUX),new=stringr::str_replace(colnames(AUX),pattern = "V","IMF-"))

IVARS <- c("t")
MVARS <- colnames(AUX[, -c("t")])
AUX2 <- melt(AUX, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
DATA.IMF <- AUX2[,.(X=t,Y=value,ID=variable)]



plot.highchart(
  color.palette ="ag_Sunset",
  yAxis.label =FALSE,
  plot.height = max(1000,100*NIMF),
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend="IMF",xAxis.legend="t",group.legend="IMF",
  yAxis.min=-OFFSET,
  data=DATA.IMF)


# ----
# Remove IMF 5,6,7 y residuos

# AUX <- data.table(X=IMF$tt,M)
# DATA[,.(t=IMF$tt,"Signal"=IMF$original.signal)]
imf_target <- c(1,2,3,4)
COLS <- paste0("V",imf_target)
DATA.IMF.R <- DT.IMF[,.(X,Y=rowSums(.SD),ID=paste0(ID_TARGET,".R")),.SDcols=COLS]
plot.ggplot2(DATA.IMF.R, plot.type = "line",line.size=0.5)


# ------
# Deteccion de modos de baja frecuencia
DT <- data.table(X=IMF$tt,M)
setnames(DT,old=colnames(DT),new=stringr::str_replace(colnames(DT),pattern = "V","IMF-"))
IVARS <- c("X")
MVARS <- colnames(DT[, -c("X")])
DT <- data.table::melt(DT, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
setnames(DT,old=c("variable","value"),new=c("ID","Y"))
FFT.IMF <- DT[,.getFFT(.SD),by="ID"]

# Plot Parameters
#
# kf <- 0.15
# fnyq <- 1/dt/2
# fmax <- round(kf*fnyq)
# DT <- FFT.IMF[f<=fmax & f>=0][,.(X=f,Y=PSD,ID=ID)]
#
# plot.highchart(
#   color.palette ="ag_Sunset",
#   plot.height = max(500,150*NIMF),
#   plot.type="spline",
#   legend.layout="horizontal",
#   legend.show=TRUE,
#   xAxis.log = FALSE,yAxis.log = FALSE,
#   yAxis.legend="PSD",xAxis.legend="Frequency[Hz]",group.legend="IMF",
#   data=DT)


# Sapply con 10 iteraciones a las tres direcciones de RAW directo\ y probar luego


