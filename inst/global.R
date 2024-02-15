
source("setup.R")



if(.Platform$OS.type=="windows"){
  RootFolder <- file.path("D:/Public/Database")
  stopifnot(dir.exists(RootFolder))
  
} 
if(.Platform$OS.type=="unix"){
  RootFolder <- file.path("~/Public/Database")
  stopifnot(dir.exists(RootFolder))
  
}
HazardFolder <- file.path(RootFolder,"gmdp","index")#"D:/Public/Database/index"
SitesFolder <- file.path(RootFolder,"dsra","index")#"D:/Public/Database/index"
SourceFolder <-  file.path(RootFolder,"gmdb","source")#"D:/Public/Database/source"
IndexFolder <-  file.path(RootFolder,"gmdb","index")#"D:/Public/Database/index"
RecordsFolder <- file.path(RootFolder,"gmdb","records")# "D:/Public/Database/records"
ExportFolder <-  file.path(RootFolder,"gmdb","export")#"D:/Public/Database/export"
NSMAX <- 200 #Sites per record
# RECOMPILE <- FALSE
# MULTI <- TRUE
# WRITE_RECORD <- FALSE
includeRawData <- FALSE
ProjectID <- "GMDB"
TargetUnits <-  "mm"
LowPassAT=TRUE
LowPassVT=TRUE
LowPassDT=TRUE
SiteResponse=TRUE
DetrendAT=TRUE
DetrendVT=TRUE
DetrendDT=TRUE
# DownFs <- 100
# GMSP <- .buildGMSP(Fs=DownFs)

FILE <- file.path(HazardFolder,"PSHA.Rds")
stopifnot(file.exists(FILE))
UserData <- readRDS(file = FILE) |> na.omit()
