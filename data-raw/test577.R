devtools::load_all()
library(xplot)

if(!exists("SET")){
  RecordsFolder <- file.path("/home/averri/Public/Database/gmdb/source/tables")
  SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
}

RSN.TARGET <- 577
ID_TARGET <- "DT.H1"

# -----
RAW <- SET[[RSN.TARGET]] #577 300 1500
RECORD <- buildTS(
  a=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  DownFs=60,
  DerivateDT=FALSE,
  DerivateVT=FALSE,
  DetrendAT=FALSE,
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

# ----
DATA <- DT.TS[ID==ID_TARGET]
# plot.highchart(DATA, plot.type = "line")
plot.ggplot2(DATA, plot.type = "line",line.size=0.5)

# -----
s <- DT.TS[ID==ID_TARGET]$Y
t <- DT.TS[ID==ID_TARGET]$X
nit <- 10
nimf <- 10
n <- -5# n<0
namp <- 0.5*10^n
ntype <- "gaussian" #c("uniform","gaussian")
IMF <- hht::CEEMD(sig=s, tt=t, noise.amp=namp, trials=nit, verbose = TRUE, noise.type=ntype)

# ----
# PLOT IMF
M <- IMF$imf
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
  plot.height = max(1000,100*nimf),
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend="IMF",xAxis.legend="t",group.legend="IMF",
  yAxis.min=-OFFSET,
  data=DATA)


# ------
# Deteccion de modos de baja frecuencia
M <- IMF$imf
DATA <- data.table(X=IMF$tt,M)
setnames(DATA,old=colnames(DT),new=stringr::str_replace(colnames(DT),pattern = "V","IMF-"))
IVARS <- c("t")
MVARS <- colnames(DT[, -c("t")])
DT <- data.table::melt(DT, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
setnames(DT,old=c("variable","value"),new=c("ID","s"))
FFTIMF <- DT[,getFFT(.SD),by="ID"]
return(FFTIMF)


# ----
# Remove IMF 5,6,7 y residuos
M <- IMF$imf
DATA <- data.table(X=IMF$tt,M)
DATA[,.(t=IMF$tt,"Signal"=IMF$original.signal)]
imf_target <- c(1,2,3)
COLS <- paste0("V",imf_target)
DATA <- DATA[,.(X,Y=rowSums(.SD),ID=paste0(ID_TARGET,".R")),.SDcols=COLS]
plot.ggplot2(DATA, plot.type = "line",line.size=0.5)



