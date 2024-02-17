devtools::load_all()

RecordsFolder <- file.path("/home/averri/Public/Database/gmdb/source/Rds")

SET <- readRDS(file.path(RecordsFolder,"AT2.Rds"))
RAW <- SET[[526]]
RECORD <- getVDA(a=RAW$AT,dt=RAW$dt,UN=RAW$SourceUnits,DownFs=100,DetrendAT=FALSE,DetrendVT=FALSE,DetrendDT=FALSE,RestoreScale=FALSE,TargetUnits="mm",NW=2048,OVLP=75)
TS <- RECORD$TS$I
