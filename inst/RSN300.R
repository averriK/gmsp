devtools::load_all()
library(data.table)
if(!exists("SET")){
  RecordsFolder <- file.path("/Users/averri/Database/gmdb/source/tables")
  SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
}

RSN_TARGET <- 300  #577 300 1500 1540
OCID_TARGET <- "UP"

# RSN 300. VT has a low frequency noise
RAW <- SET[[RSN_TARGET]] #577 300 1500

# -----
# Stage 1/ Raw record

R1 <- buildTS(
  x=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  Order=2,
  Fmax=15,
  Resample = FALSE,
  LowPass = TRUE,
  TargetUnits="mm",
  removeIMF1 = 0,
  removeIMFn = 0)
TSL <- R1$TSL

DATA <- TSL[OCID==OCID_TARGET,.(X=t,Y=s,ID=ID)]
xplot::plot.highchart(
  color.palette ="Dynamic",
  yAxis.label =TRUE,
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  xAxis.legend="t",
  data=DATA)

# -----
# Stage 2/ Remove Residues and 1 IMF (low freq) from AT
devtools::load_all()

R2 <- buildTS(
  x=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  Order=2,
  Fmax=15,
  Resample = TRUE,
  LowPass = TRUE,
  TargetUnits="mm",
  removeIMF1 = 0,
  removeIMFn = 4)
TSL <- R2$TSL

DATA <- TSL[OCID==OCID_TARGET,.(X=t,Y=s,ID=ID)]
xplot::plot.highchart(
  color.palette ="Dynamic",
  yAxis.label =TRUE,
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  xAxis.legend="t",
  data=DATA)
# VT records takes too long... it seems that it takes too long to remove the 1st IMF
