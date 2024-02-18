rm(list=ls())
devtools::load_all()
library(xplot)
RecordsFolder <- file.path("/home/averri/Public/Database/gmdb/source/tables")
SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
RAW <- SET[[1000]]
RECORD <- getVDA(
  a=RAW$AT,
  dt=RAW$dt,
  UN=RAW$SourceUnits,
  DownFs=200,
  DerivateDT=FALSE,
  DerivateVT=FALSE,
  DetrendAT=TRUE,
  DetrendVT=FALSE,
  DetrendDT=TRUE,
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
# plot.highchart(DT.TS[ID=="AT.H1"], plot.type = "line")
# plot.highchart(DT.TS[ID=="VT.H1"], plot.type = "line")
# plot.highchart(DT.TS[ID=="DT.H1"], plot.type = "line")


s <- DT.TS[ID=="DT.H1"]$Y
t <- DT.TS[ID=="DT.H1"]$X
nit <- 20
nimf <- 10
n <- -7# n<0
namp <- 0.5*10^n
ntype <- "gaussian" #c("uniform","gaussian")
IMF <- hht::CEEMD(sig=s, tt=t, noise.amp=namp, trials=nit, verbose = TRUE, noise.type=ntype)

M <- IMF$imf
OFFSET <- 1.25*ceiling(max(M)-min(M))
for(i in 1:ncol(M)){
  j <- ncol(M)-i+1
  M[,j] <- M[,j]+OFFSET*i
}
DT.IMF <- as.data.table(M)
DT.IMF[,`:=`(t=IMF$tt,"Residue"=IMF$residue,"Signal"=IMF$original.signal+OFFSET*(ncol(M)+2))]
setnames(DT.IMF,old=colnames(DT.IMF),new=stringr::str_replace(colnames(DT.IMF),pattern = "V","IMF-"))

IVARS <- c("t")
MVARS <- colnames(DT.IMF[, -c("t")])
DT.IMF <- melt(DT.IMF, id.vars = IVARS, measure.vars = MVARS) |> na.omit()
DT.IMF <- DT.IMF[,.(X=t,Y=value,ID=variable)]
plot.highchart(
  color.palette ="ag_Sunset",
  yAxis.label =FALSE,
  plot.height = max(1000,100*nimf),
  plot.type="line",
  legend.layout="horizontal",
  legend.show=TRUE,
  yAxis.legend="IMF",xAxis.legend="t",group.legend="IMF",
  yAxis.min=-OFFSET,
  data=DT.IMF)





