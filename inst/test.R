devtools::load_all()

if(!exists("SET")){
  RecordsFolder <- file.path("/Users/averri/Database/gmdb/source/tables")
  SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
}
RSN_TARGET <- 300  #577 300 1500 1540
# RSN 300. VT has a low frequency noise
RAW <- SET[[RSN_TARGET]] #577 300 1500

# -----
# Stage 1/ Raw record

R1 <- buildTS(
  x=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  Fmax=25,
  TrimZeros = TRUE,
  Rebuild = FALSE,
    Resample = TRUE,
  LowPass.AT = TRUE,
  LowPass.VT = FALSE,
  LowPass.DT = FALSE,
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
ID_TARGET <- "VT"
OCID_TARGET <- "UP"
DATA <- TSL[ID==ID_TARGET & OCID==OCID_TARGET,.(X=t,Y=s,ID=paste0(ID,".",OCID))]
xplot::plot.highchart(
  color.palette ="Blue-Red",
  yAxis.label =TRUE,
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend=paste0(ID_TARGET,".",OCID_TARGET),xAxis.legend="t",
  data=DATA)

# check imfs
TS <- TSL[ID==ID_TARGET & OCID==OCID_TARGET]
AUX <- buildIMF(t=TS$t,s=TS$s,method="emd",plot=TRUE)
DATA <- AUX$plot.data
offset <- AUX$plot.offset
nimf <- AUX$nimf
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

# -----
# Stage 2/ Remove Residues and 1 IMF (low freq) from AT

R2 <- buildTS(
  x=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  Fmax=15,
  LowPass.AT = TRUE,
  LowPass.VT = TRUE,
  LowPass.DT = TRUE,
  Resample = TRUE,
  TrimZeros = TRUE,
  Rebuild = TRUE,
  TargetUnits="mm",
  EMD.method="emd",
  EMD.AT = TRUE,
  EMD.VT = FALSE,
  EMD.DT = TRUE,
  removeIMF1.AT = 0,# 0
  removeIMFn.AT = 1,# 1:
  removeIMF1.VT = 0,
  removeIMFn.VT = 0,
  removeIMF1.DT = 0,#
  removeIMFn.DT = 1, # 1: works
  NW=2048,
  OVLP=75)
TSL <- R2$TSL
ID_TARGET <- "DT"
OCID_TARGET <- "UP"
DATA <- TSL[ID==ID_TARGET & OCID==OCID_TARGET,.(X=t,Y=s,ID=paste0(ID,".",OCID))]
xplot::plot.highchart(
  color.palette ="Blue-Red",
  yAxis.label =TRUE,
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend=paste0(ID_TARGET,".",OCID_TARGET),xAxis.legend="t",
  data=DATA)
# VT records takes too long... it seems that it takes too long to remove the 1st IMF
