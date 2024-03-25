devtools::load_all()

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
  Taper=3,
  Fmax=15,
  TrimZeros = TRUE,
  Resample = FALSE,
  Detrend.AT = TRUE,
  Detrend.VT = TRUE,
  Detrend.DT = TRUE,
  LowPass.AT = TRUE,
  LowPass.VT = TRUE,
  LowPass.DT = TRUE,
  TargetUnits="mm",
  EMD.method="emd",
  EMD.AT = FALSE,
  EMD.VT = FALSE,
  EMD.DT = FALSE,
  removeIMF1.AT = 0,
  removeIMFn.AT = 0,
  removeIMF1.VT = 0,
  removeIMFn.VT = 0,
  removeIMF1.DT = 0,
  removeIMFn.DT = 0,
  NW=2048,
  OVLP=75)
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
  Fmax=15,
  Taper=3,
  Detrend.AT = TRUE,
  Detrend.VT = TRUE,
  Detrend.DT = TRUE,
  LowPass.AT = TRUE,
  LowPass.VT = TRUE,
  LowPass.DT = TRUE,
  Rebuild = TRUE,
  Resample = FALSE,
  TrimZeros = TRUE,
  TargetUnits="mm",
  EMD.method="emd",
  EMD.AT = TRUE,
  EMD.VT = TRUE,
  EMD.DT = TRUE,
  removeIMF1.AT = 1,# 0
  removeIMFn.AT = 0,# 1:
  removeIMF1.VT = 0,
  removeIMFn.VT = 1,
  removeIMF1.DT = 0,#
  removeIMFn.DT = 1, # 1: works
  NW=2048,
  OVLP=75)
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
