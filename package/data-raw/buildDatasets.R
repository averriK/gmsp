library(data.table)
# Scale A in g, V in cm/s
go <- 9.88 #[m/s2]
mm2m2g <- go*1000 #[mm/m]
mm2cm <- 10
mm2m <- 1000

PATH <-  "/home/averri/Public/Database/gmdb/index"
DT <- readRDS(file=file.path(PATH,"OwnerIntensityTable.Rds"))
DT <- DT[PGA.cov<=0.1 & IA.cov<=0.1 & ARMS.cov<=0.1]
DT <- DT[PGA/ARMS<7]
# DT <- DT[(isUP==TRUE & PGA>50 & PGA < 10000)| isUP==FALSE]
DT <- DT[(isUP==FALSE & PGA>50 & PGA < 15000)| isUP==TRUE]

OIT <- DT[UN=="mm" & TopLayerID=="station",.(
  RecordSN,
  D0595,Dmax,TmA,
  PGAs=PGA/mm2m2g,PGVs=PGV/mm2m,PGDs=PGD/mm2m,
  ARMSs=ARMS/mm2m2g,VRMSs=VRMS/mm2m,DRMSs=DRMS/mm2m,
  CAV5s=CAV5/mm2m,IAs=IA/mm2m,
  M=Magnitude,
  Re=EpicenterDistance,
  S=ifelse(OwnerNEHRP %in% c("A","B","BC"),1,ifelse(OwnerNEHRP %in% c("C","CD"),2,ifelse(OwnerNEHRP %in% c("D","DE","E"),3,0))) |> as.factor(),
  V=isUP,EpicenterLongitude,EpicenterLatitude)][] |> na.omit()


RecordsFolder <- file.path("/home/averri/Public/Database/gmdb/source/Rds")
DIRS <- list.dirs(path = RecordsFolder, full.names = FALSE, recursive = FALSE)
# DIR <- DIRS[3]
# DIRS <- c("V2A","TRZ","TRA","TRB")
for(DIR in DIRS){
  # Load Data
  IPATH <- file.path(RecordsFolder,DIR)
  RECORDS <- list.files(path = IPATH, pattern = "Rds", full.names = FALSE) |> stringr::str_remove(pattern=".Rds")
  RECORDS <- RECORDS[RECORDS %in% OIT$RecordSN]
  if(length(RECORDS)==0) next()
  FILES <- file.path(IPATH,paste0(RECORDS,".Rds"))
  DATA <- list()
  # FILE <- FILES[11]
  for(FILE in FILES){
    # Load Data
    RAW <- readRDS(FILE)
    AT <- RAW$AT
    ID <- RAW$RecordSN
    OCID <- colnames(AT)
    UP_ID <- c("UD","UP","V","HLZ","HNZ","DWN","DN","HHZ","Z","DOWN","VER","VERT","VERTICAL","VRT","BHZ","HHZ","HNZ","HLZ","HNZ","HLZ","HN3","HGZ","UD","UP","V","UPDO") |> unique()
    idx <- toupper(OCID)%in%UP_ID
    if(sum(idx)!=1){
      next() # invalid record. No vertical records or multiple matches for UP_ID
    }
    setnames(AT,old=OCID[idx],new="UP")
    COLS <- colnames(AT[,-"UP"])

    PGA <- AT[,sapply(.SD,function(x){max(abs(x))}),.SDcols = -c("UP")]
    idx <- PGA == max(PGA)
    if(sum(idx)==2) {
      setnames(AT,old=COLS[1],new="H1")
      setnames(AT,old=COLS[2],new="H2")
    } else if(sum(idx)==1) {
      setnames(AT,old=COLS[idx],new="H1")
      setnames(AT,old=COLS[!idx],new="H2")
    } else {
      next()}
    RECORD <- list()
    RECORD[[ID]] <- list(RecordSN=ID,Units=RAW$Units,dt=RAW$dt,Header=RAW$Header,OwnerRecordSN=RAW$OwnerRecordSN,AT=AT,OCID=OCID)

    DATA <- c(DATA, RECORD)

  }
  # Save Data
  saveRDS(DATA, file = file.path("data-raw",paste0(DIR,".Rds")))
}
