devtools::load_all()

RecordsFolder <- file.path("/home/averri/Public/Database/gmdb/source/Rds")

SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
RAW <- SET[[526]]
RECORD <- getVDA(a=RAW$AT,dt=RAW$dt,UN=RAW$SourceUnits,DownFs=100,DetrendAT=FALSE,DetrendVT=FALSE,DetrendDT=FALSE,RestoreScale=FALSE,TargetUnits="mm",NW=2048,OVLP=75)
TS <- RECORD$TS$I
dt <- RECORD$dt
ts <- seq(0,dt*(nrow(TS)-1),dt)
TS <- data.table(ts=ts,TS)
# saveRDS((RAW$AT),file.path(RecordsFolder,"data-raw/AT.Rds"))

library(xplot)
DATA <- rbindlist(
  list(
    TS[,.(X=ts,Y=TS$AT.H1,ID="AT")],
    TS[,.(X=ts,Y=TS$VT.H1,ID="VT")],
    TS[,.(X=ts,Y=TS$DT.H1,ID="DT")]
  )
)

plot.highchart(DATA, plot.type = "line")
