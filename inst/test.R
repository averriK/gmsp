devtools::load_all()
library(xplot)
library(data.table)
if(!exists("SET")){
  RecordsFolder <- file.path("/Users/averri/Database/gmdb/source/tables")
  SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))

}

RSN.TARGET <- 1500  #577 300 1500 1540
ID_TARGET <- "AT.H1"
COMPLETE <- FALSE

# -----
RAW <- SET[[RSN.TARGET]] #577 300 1500
RECORD <- buildTS(
  a=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  DownFs=250,
  FlatZerosAT=FALSE,
  DerivateDT=FALSE,
  DerivateVT=FALSE,
  DetrendAT=TRUE,
  DetrendVT=FALSE,
  DetrendDT=FALSE,
  RestoreScale=FALSE,
  TargetUnits="mm",
  NW=2048,
  OVLP=75)
TS <- RECORD$TS$I
dt <- RECORD$dt
ts <- seq(0,dt*(nrow(TS)-1),dt)
TS <- data.table(ts=ts,TS)
# saveRDS((RAW$AT),file.path(RecordsFolder,"data-raw/AT.Rds"))
IVARS <- c("ts","W")
MVARS <- colnames(TS[, -c("ts","W")])
DT.TS <- data.table::melt(TS, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
DT.TS <- DT.TS[,.(X=ts,Y=value,ID=variable)]
DATA <- DT.TS[ID==ID_TARGET]
# plot.highchart(DATA, plot.type = "line")
plot.ggplot2(DATA, plot.type = "line",line.size=0.5)

# -----
s <- DT.TS[ID==ID_TARGET]$Y
t <- DT.TS[ID==ID_TARGET]$X
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

## Cheaper version


# ----
# PLOT IMF
OFFSET <- 1.25*ceiling(max(M)-min(M))
for(i in 1:ncol(M)){
  j <- ncol(M)-i+1
  M[,j] <- M[,j]+OFFSET*i
}
DATA <- as.data.table(M)
DATA[,`:=`(t=IMF$tt,"Residue"=IMF$residue,"Signal"=IMF$original.signal+OFFSET*(ncol(M)+2))]
setnames(DATA,old=colnames(DATA),new=stringr::str_replace(colnames(DATA),pattern = "V","IMF-"))

IVARS <- c("t")
MVARS <- colnames(DATA[, -c("t")])
DATA <- melt(DATA, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
DATA <- DATA[,.(X=t,Y=value,ID=variable)]
plot.highchart(
  color.palette ="ag_Sunset",
  yAxis.label =FALSE,
  plot.height = max(1000,100*NIMF),
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend="IMF",xAxis.legend="t",group.legend="IMF",
  yAxis.min=-OFFSET,
  data=DATA)


# ----
# Remove IMF 5,6,7 y residuos

DATA <- data.table(X=IMF$tt,M)
DATA[,.(t=IMF$tt,"Signal"=IMF$original.signal)]
imf_target <- c(2,3)
COLS <- paste0("V",imf_target)
DATA <- DATA[,.(X,Y=rowSums(.SD),ID=paste0(ID_TARGET,".R")),.SDcols=COLS]
plot.ggplot2(DATA, plot.type = "line",line.size=0.5)


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


