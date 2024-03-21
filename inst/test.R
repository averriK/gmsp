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
  TrimZeros=TRUE,
  DetrendAT=TRUE,
  DetrendVT=TRUE,
  DetrendDT=TRUE,
  PadZeros=TRUE,
  TargetUnits="mm",
  NW=1024,
  OVLP=75)
TSL <- RECORD$TSL
TSW <- RECORD$TSW
dt <- RECORD$dt


DATA <- TSL[ID=="AT" & OCID=="H1",.(X=t,Y=s,ID=paste0(ID,".",OCID))]
plot.ggplot2(DATA, plot.type = "line",line.size=0.5)

# -----
s <- DATA$Y
t <- DATA$X
IMF <- .getEMD(t=t,s=s)
PLOT <- .plotEMD(IMF)
PLOT
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


