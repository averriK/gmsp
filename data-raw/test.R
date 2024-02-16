AT2 <- readRDS("data-raw/AT2.Rds")[[11]]
dt <- AT2$dt
AT <- AT2$AT
SourceUnits <- AT2$Units
getVDA(a=AT,dt=dt,UN=SourceUnits,DownFs=100,DetrendAT=FALSE,DetrendVT=FALSE,DetrendDT=FALS,RestoreScale=FALSE,TargetUnits="mm",NW=2048,OVLP=75)
