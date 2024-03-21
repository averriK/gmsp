devtools::load_all()
library(xplot)
library(data.table)
if(!exists("SET")){
  RecordsFolder <- file.path("/Users/averri/Database/gmdb/source/tables")
  SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
}

RSN_TARGET <- 577  #577 300 1500 1540
ID_TARGET <- "DT"
OCID_TARGET <- "H1"
COMPLETE <- FALSE

# -----
RAW <- SET[[RSN_TARGET]] #577 300 1500
RECORD <- buildTS(
  x=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  Fmax=20,
  RebuildAT = TRUE,
  TargetUnits="mm",
  RemoveFirstIMF = TRUE,
  RemoveLastIMF = TRUE,
  NW=1024,
  OVLP=75)
TSL <- RECORD$TSL
TSW <- RECORD$TSW
dt <- RECORD$dt


# DT.TS <- TSL[ID==ID_TARGET & OCID==OCID_TARGET,.(X=t,Y=s,ID=paste0(ID,".",OCID))]
DT.TS <- TSW[,.(X=ts,Y=AT.H1,ID="AT.H1")]
plot.ggplot2(DT.TS, plot.type = "line",line.size=0.5)
s <- DT.TS$Y
t <- DT.TS$X
# -----

AUX <- .buildIMF(t=t,s=s)
IMF <- AUX$imf
nimf <- ncol(IMF)
sR <- IMF[,-c(1),with = FALSE][,rowSums(.SD)]
DATA <- rbindlist(list(DT.TS,data.table(X=t,Y=sR,ID="Filtered")))
xplot::plot.highchart(
  color.palette ="Blue-Red",
  yAxis.label =TRUE,
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend=paste0(ID_TARGET,".",OCID_TARGET),xAxis.legend="t",
  data=DATA)


DATA <- AUX$plot.data
offset <- AUX$plot.offset
xplot::plot.highchart(
  color.palette ="ag_Sunset",
  yAxis.label =FALSE,
  plot.height = max(1000,100*nimf),
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend="IMF",xAxis.legend="t",group.legend="IMF",
  yAxis.min=offset,
  data=DATA)
# ----


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


