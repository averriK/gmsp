# rm(list=ls())
library(readr)
library(data.table)
library(stringr)
library(triangle)
library(digest)
library(seewave)
library(pracma)
library(purrr)
library(lubridate)
library(anytime)
library(shiny)
library(ggplot2)
library(GGally)
library(ggpubr)
library(compiler)#library(compiler, lib.loc = "/usr/lib64/R/library")
options(readr.show_progress = FALSE)
options(readr.num_columns = 0)
options(rstudio.help.showDataPreview = FALSE)


## -- 1. Ground-Motion File Parser (GMFP)  -------------------------------------------------
# Copyright SRK (C) 2010-2021
# Author: A.Verri K. (averri@srk.com.ar)
# Maintainer: A.Verri K. (averri@srk.com.ar)
# Contributors: M. Balbi,  J. Mussatt
# Date: 26-05-2021

.setParseOptions <- function(TypeID=NA){
  stopifnot(!is.na(TypeID))
  ParseOptions <- switch(
    TypeID,
    ## -----------------------------------------------------
    "ACA"=list(
      TYPE=TypeID,
      EXT = "*.aca",
      HEADER_LINES=NA,
      COLUMN_NAMES=TRUE,
      NCOLS=NA,
      PATTERN_BLOCK_START = NA,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR=NA,
      # PATTERN_DT ="(?<=(?i)MUESTREO).*(?=(?i)m)",# remover \\:
      PATTERN_DT ="(?<=(?i)MUESTREO).*(?=(?i)m)",# remover \\:
      PATTERN_NP ="(?<=(?i)MUESTRAS).*",# remover \\:
      PATTERN_UN ="(?<=(?i)UNIDADES).*(?=/(?i)s)"#remover \\:
    ),
    ## -----------------------------------------------------
    "ACB"=list(
      TYPE=TypeID,EXT = "*.acb",
      HEADER_LINES=NA,
      COLUMN_NAMES=TRUE,
      NCOLS=NA,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR=NA,
      PATTERN_DT ="(?<=\\((?i)Hz\\)).*",# remover \\:
      PATTERN_NP ="(?<=(?i)SAMPLES).*",# remover \\:
      PATTERN_UN ="(?<=(?i)UNITS).*(?=/(?i)s)"#remover \\:
    ),
    ## -----------------------------------------------------
    "LIS"=list(
      TYPE=TypeID,
      EXT = "*.lis",
      HEADER_LINES=NA,
      COLUMN_NAMES=TRUE,
      NCOLS=NA,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR=NA,
      PATTERN_DT ="(?<=(?i)Delta\\s(?i)t).*",# remover \\:
      PATTERN_NP ="(?<=(?i)points).*",# remover \\:
      PATTERN_UN ="(?<=(?i)DATA\\s(?i)IN).*(?[[:punct:]]+)"#remover \\:
    ),
    ## -----------------------------------------------------
    "TRB"=list(
      TYPE=TypeID,
      EXT = "*.trb",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=1,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR="(?<=(?i)Componente\\:).*",# remover \\s
      PATTERN_DT ="(?<=(?i)muestreo\\:).*(?=(?i)muestras)",# remover \\s
      PATTERN_NP ="(?<=(?i)muestras\\:).*",# remover \\s
      PATTERN_UN ="(?<=(?i)Unidades\\:).*(?=/(?i)s)"
    ),
    ## -----------------------------------------------------
    "TRC"=list(
      TYPE=TypeID,
      EXT = "*.trc",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=1,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR="(?<=(?i)STREAM\\:).*",# remover \\s
      PATTERN_DT ="(?<=(?i)SAMPLING_INTERVAL_S\\:).*",# remover \\s
      PATTERN_NP ="(?<=(?i)NDATA\\:).*",# remover \\s
      PATTERN_UN ="(?<=(?i)UNITS\\:).*(?=/(?i)s)"
    ),
    ## -----------------------------------------------------
    "TRA"=list(
      TYPE=TypeID,
      EXT = "*.tra",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=15,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR="(?<=(?i)Component\\:).*",# remover \\s
      PATTERN_DT ="(?<=(?i)rate\\:).*",# remover \\s
      PATTERN_NP ="(?<=(?i)padding\\:).*",# remover \\s
      PATTERN_UN =".*(?=/(?i)s)"
    ),
    ## -----------------------------------------------------
    "V2A"=list(
      TYPE=TypeID,EXT = "*.v2a",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=10,
      PATTERN_BLOCK_START = "^(?i)Corrected",
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR= "(?<=(?i)Component\\s)[[:alnum:]]+" , #IDX_DIR= 13 ,
      PATTERN_DT ="(?<=(?i)corrected\\sdata\\sat\\s).*(?=\\s(?i)sec)",#IDX_DT=11,
      PATTERN_NP ="(?<=(?i)Number\\sof\\spoints\\s).*(?=\\s*(?i)D)",#IDX_NP=10,
      PATTERN_UN ="[[:alpha:]]+(?=/(?i)s/(?i)s)"
    ),
    ## -----------------------------------------------------
    "V2C"=list(
      TYPE=TypeID,
      EXT = "*.v2c",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=8,
      PATTERN_BLOCK_START = "^(?i)CORRECTED",
      PATTERN_TS_START="\\.*ACCEL\\sDATA",
      PATTERN_TS_END ="/&",
      PATTERN_DIR= "(?<=\\:).*(?=\\((?i)STA)" ,
      PATTERN_DT ="(?<=(?i)INTERVALS\\sOF).*(?=(?i)S)",
      PATTERN_NP =".*(?=POINTS\\sOF\\sINSTRUMENT)",
      PATTERN_UN ="[[:alpha:]]+(?=/(?i)SEC/(?i)SEC)"
    ),
    ## -----------------------------------------------------
    "AT2"=list(
      TYPE=TypeID,EXT = "*.at2",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=8,
      PATTERN_BLOCK_START = NA,
      PATTERN_TS_START= NA,
      PATTERN_TS_END = NA,
      PATTERN_DIR= NA , #split ,
      PATTERN_DT = "(?<=(?i)DT\\=).*(?=(?i)S)",
      PATTERN_NP ="(?<=(?i)NPTS\\=)\\s*[[:digit:]]+",
      PATTERN_UN ="(?<=(?i)UNITS\\sOF).*"#IDX_UN=46
    ),
    ## -----------------------------------------------------
    "TRZ"=list(
      TYPE=TypeID,
      EXT = "*.TRZ",
      HEADER_LINES=NA,
      COLUMN_NAMES=FALSE,
      NCOLS=15,
      PATTERN_TS_START=NA,
      PATTERN_TS_END =NA,
      PATTERN_DIR="(?<=(?i)Component\\:).*",# remover \\s
      PATTERN_DT ="(?<=(?i)rate\\:).*",# remover \\s
      PATTERN_NP ="(?<=(?i)padding\\:).*",# remover \\s
      PATTERN_UN =".*(?=/(?i)s)"
    )
  )
  return(ParseOptions)

}

.buildFileTable <- function(TypeID,MULTI=FALSE,SourceFolder){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  stopifnot(toupper(TypeID) %in% c("TRZ","TRA","TRB","TRC","V2A","ACA","ACB","LIS","V2C","AT2"))
  # Build FileTable - ------------------------------------------------------
  # Set Input/Output Paths
  PATH <- file.path(SourceFolder,"raw",TypeID)
  # Get EXT
  EXT <- .setParseOptions(TypeID)$EXT

  cat(sprintf("> Building %s Index from %s...",TypeID,PATH))
  # Build Filetable
  FT <- data.table(
    OwnerFilename=list.files(
      path = PATH,
      pattern = EXT,
      recursive = TRUE,
      ignore.case=TRUE,
      include.dirs = FALSE)
  )

  # include relative path
  FT[,OwnerFilename:=file.path(PATH,OwnerFilename)]

  # Remove EXT
  FT[,OwnerFileSN:=substr(basename(OwnerFilename),start = 1,stop = nchar(basename(OwnerFilename))-4)]

  #  Remove Duplicates
  idx <- duplicated(FT$OwnerFileSN)
  if(sum(idx)>0) FT <- FT[!idx]
  FT[,TypeID:=TypeID]
  FT[,Duplicated:=duplicated(OwnerFileSN)]

  # Filename does not include ComponentID
  # OwnerRecordSN is OwnerFileSN
  if(TypeID %in% c("ACA","ACB","V2C","V2A","LIS")){
    FT[,is3D:=TRUE]
    FT[,OwnerComponentID:="*"]


    if(TypeID=="ACA"){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[_]",simplify = TRUE)  |>
        as.data.table()
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V1,AUX$V2,AUX$V3) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerStationSN:=AUX$V3 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:=paste0(AUX$V1,AUX$V2) |>
      #      gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]

      # FT[,OwnerRecordSN:=paste0(OwnerEventSN,OwnerStationSN) ] # Remove Station ID
    }

    if(TypeID=="ACB"){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[_]",simplify = TRUE)  |>
        as.data.table()
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V2,AUX$V1) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerStationSN:=AUX$V1 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:=AUX$V2 |>
      #      gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
    }

    if(TypeID=="V2A"){
      FT[,OwnerRecordSN:=paste0(TypeID,basename(OwnerFileSN))|> gsub(pattern = "[[:punct:]]", replacement = "")]
      # FT[,OwnerEventSN:=NA]
      FT[,OwnerStationSN:=NA]
    }

    if(TypeID=="V2C"){
      FT[,OwnerRecordSN:=paste0(TypeID,basename(OwnerFileSN))|> gsub(pattern = "[[:punct:]]", replacement = "")]
      # FT[,OwnerEventSN:=NA]
      FT[,OwnerStationSN:=NA]
    }

    if(TypeID=="LIS"){
      AUX <- data.table(V1=FT$OwnerFileSN)
      AUX[,V2:=substr(V1,nchar(V1)-3,nchar(V1))]
      AUX[,V3:=paste0(substr(V1,1,nchar(V1)-4),"00")]
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V3,AUX$V2) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerStationSN:=AUX$V2 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:=AUX$V3 |>
      #      gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
    }
  }

  # Filename includes ComponentID
  # OwnerRecordSN is OwnerFileSN without OwnerComponentID
  if(TypeID %in% c("AT2","TRA","TRB","TRC","TRZ")){
    FT[,is3D:=FALSE]
    FT[,OwnerComponentID:=NA]

    if(TypeID=="AT2"){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[_]",simplify = TRUE)  |>
        as.data.table()
      AUX[,V4:=substr(V3,nchar(V3)-2,nchar(V3))]
      AUX[,V5:=substr(V3,1,nchar(V3)-3)]
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V1,AUX$V2,AUX$V5) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerStationSN:=AUX$V5 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:=NA] #paste0(AUX$V1,AUX$V2,AUX$V5,AUX$V6)
      FT[,OwnerComponentID:=  AUX$V4  |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]

    }

    if(TypeID %in% c("TRA")){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[.]",simplify = TRUE)  |>
        as.data.table()
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V1,AUX$V2,AUX$V3,AUX$V4,AUX$V5,AUX$V6,AUX$V7) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:= NA]
      FT[,OwnerStationSN:=AUX$V7|>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerComponentID:= AUX$V8 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
    }


    if(TypeID %in% c("TRZ")){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[.]",simplify = TRUE)  |>
        as.data.table()
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V1,AUX$V2,AUX$V3,AUX$V4,AUX$V5,AUX$V6,AUX$V7,AUX$V8,AUX$V10) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:= paste0(AUX$V1,AUX$V2,AUX$V3,AUX$V4,AUX$V5,AUX$V6,AUX$V7)]
      FT[,OwnerStationSN:=AUX$V8|>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerComponentID:= AUX$V9 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
    }

    if(TypeID=="TRB"){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[-]",simplify = TRUE) |> as.data.table()
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V1,AUX$V2,AUX$V3) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:=paste0(AUX$V1,AUX$V2) |>
      #      gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerStationSN:=AUX$V3 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerComponentID:= AUX$V4 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]

    }

    if(TypeID=="TRC"){
      AUX <- str_split(string=FT$OwnerFileSN,pattern="[.]",simplify = TRUE) |> as.data.table()
      FT[,OwnerRecordSN:=paste0(TypeID,AUX$V1,AUX$V2,AUX$V4,AUX$V5,AUX$V6,AUX$V7,AUX$V8) |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      # FT[,OwnerEventSN:=paste0(AUX$V5,AUX$V6) |>
      #      gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerStationSN:=AUX$V2 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
      FT[,OwnerComponentID:=AUX$V3 |>
           gsub(pattern="[[:punct:]]",replacement = "") |> toupper()]
    }
  }

  cat(sprintf("Done!\n"))

   # Update FileTable - ------------------------------------------------------
  PATH <- file.path(SourceFolder,"Rds",TypeID)
  cat(sprintf("> Updating %s Index from %s...",TypeID,PATH))
  # Set Serial
  FT[,RecordSN:=sapply(basename(OwnerRecordSN),function(x){toupper(digest(object=x,algo="crc32", serialize=FALSE))})]
  FT[,SourceFilename:=file.path(PATH,paste0(RecordSN,".Rds"))]

  # FT[,BlockSN:=str_sub(RecordSN,start = 1L, end = 2L)]
  # FT[,RecordFilename := file.path(RecordsFolder,"Rds",BlockSN,basename(SourceFilename))]

  # FT[, Regularized := RecordFilename |> file.exists()]
  FT[, Processed := SourceFilename |> file.exists()]

  # MULTI (useful for multiple session)
  if(MULTI==TRUE){FT <- FT[sample(1:nrow(FT))]}
  # else {setorder(FT,RecordSN)}

  cat(sprintf("Done!\n"))
  # Save --------------------------------------------------------
  # cat(sprintf("> Writting File Table ..."))
  # PATH <- file.path(IndexFolder)
  # RDSfile <- file.path(PATH,paste0(TypeID,".FileTable.Rds"))
  # saveRDS(file=RDSfile,compress=TRUE, object=FT)
  # cat(sprintf("Done!\n"))
  return(FT[])
}

.parseOwnerRecordsDatabase <- function(TypeID,MULTI=FALSE, GMSP=NULL, OVERRIDE=FALSE){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  stopifnot(toupper(TypeID) %in% c("TRZ","TRA","TRB","TRC","V2A","ACA","ACB","LIS","V2C","AT2") && !is.null(GMSP))
  # Set Input/Output Paths
  RPATH <- file.path(SourceFolder,"raw",TypeID)
  # stopifnot(dir.exists(RPATH))

  OPATH <- file.path(SourceFolder,"Rds",TypeID)
  if(!dir.exists(OPATH)){dir.create(OPATH,showWarnings = FALSE)}

  # Build FileTable from Raw Files
  FT <- .buildFileTable(TypeID=TypeID,MULTI=MULTI)

  # Subset FileTable
  if(OVERRIDE==FALSE){
    SFT <- FT[Processed == FALSE & Duplicated == FALSE]
  } else {
    SFT <- FT[ Duplicated == FALSE]
  }

  RSN <- SFT$RecordSN |> unique()
  # NF <- nrow(SFT)
  NF <- length(RSN)
  for(i in seq_len(length.out = NF)){
    SET <- SFT[RecordSN == RSN[i]]
    OK <- (nrow(SET)==1 & SET$is3D[1]==TRUE) | SET$is3D[1]==FALSE
    stopifnot(OK)
    cat(sprintf("> Parse record #%s in database %s (%d/%d):",RSN[i],TypeID,i,NF))
    SourceFilename <- SET$SourceFilename |> unique()
    if(!file.exists(SourceFilename)|| OVERRIDE == TRUE){
      Record <- NULL
      # 3D RECORDS
      if(TypeID %in% c("ACA","ACB","LIS") & nrow(SET)==1 ){
        Record <- .readAC(SET)
      }
      if(TypeID %in% c("V2A","V2C") & nrow(SET)==1 ){
        Record <- .readV2(SET)
      }
      # 1D RECORDS
      if(TypeID=="AT2" & nrow(SET)==3 ){
        Record <- .readAT(SET)
      }
      if(TypeID %in% c("TRZ","TRA","TRB","TRC") & nrow(SET)==3 ) {
        Record <- .readTR(SET)
      }
      if(!is.null(Record)){
        saveRDS(compress=FALSE, object=Record,file=SourceFilename)
        cat(sprintf(" Done!.\n"))
      } else {
        cat(sprintf(" * Broken Record. Skipped.\n"))}
    } else {
      cat(sprintf("( * Processed already. Skipped.\n"))
    }
  }


}

.readTR <- function(FT=data.table()){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  rawGMSR <- NULL
  OK <- is.data.table(FT) && FT$TypeID[1] %in% c("TRZ","TRA","TRB","TRC")
  stopifnot(OK)
  # tryCatch(
  #   expr = {

      DIR <- vector(mode = "character",length = nrow(FT))
      dt <- vector(mode = "double",length = nrow(FT))
      NP <- vector(mode = "double",length = nrow(FT))
      UN <- vector(mode = "character",length = nrow(FT))
      TSL <- list()

      for(k in seq_len(length.out = nrow(FT))){
        cat(sprintf(" %s ",FT$OwnerFileSN[k]))
        TypeID <- FT$TypeID[k]
        OwnerRecordSN <- FT$OwnerRecordSN[k]
        OwnerStationSN <- FT$OwnerStationSN[k]
        # OwnerEventSN <- FT$OwnerEventSN[k]
        RecordSN <- FT$RecordSN[k]
        ParseOptions <- .setParseOptions(TypeID)
        stopifnot(!is.null(ParseOptions))
        OwnerFilename <- FT$OwnerFilename[k]

        RAW <- read_lines(OwnerFilename,skip_empty_rows = FALSE) |> str_squish()
        RAW <- RAW[RAW!=""]
        HL <- NULL
        if(TypeID %in% c("TRZ","TRA")){
          HL <- grep(x=RAW,pattern = "END_HEADER")
        }
        if(TypeID == "TRB"){
          HL <- RAW |> grep(pattern = "Unidades")
        }

        if(TypeID == "TRC"){
          HL <- RAW |> grep(pattern = "USER5")
        }
        # stopifnot(length(HL)==1 & HL>0)
        if(length(HL)!=1 | HL==0){return(NULL)}

        HEADER <- RAW[c(1:HL)]
        BLOCK <- fread(text=RAW[HL+3:length(RAW)],fill = TRUE) #RAW[-c(1:HL)]
        # Remove rows with all NAs.
        # DO NOT USE  na.omit() since could be a full column with NAs
        BLOCK <- BLOCK[rowSums(is.na(BLOCK))<2,]
        AUX <- str_match_all(string =HEADER,pattern =  ParseOptions$PATTERN_DIR) |>
          unlist()  |>  str_remove_all(pattern = "\\s*")
        DIR[k] <- AUX[AUX!=""]
        if( length(DIR[k])!=1){return(NULL)}

        dt[k] <- str_match_all(string =HEADER,pattern =  ParseOptions$PATTERN_DT) |>
          unlist()  |>  str_remove_all(pattern = "\\s*") |> as.numeric()
        # stopifnot( length(dt[k])==1 & dt[k]>0)
        if(length(dt[k])!=1 | dt[k]==0){return(NULL)}
        if(TypeID %in% c("TRA","TRZ","TRB")){
          dt[k] <- 1/dt[k]
        }

        if(TypeID %in% c("TRA","TRZ")){
          # UN_TEXT <- as.vector(RAW[HL+2] |> str_split(pattern = " ",simplify = TRUE)) |> last()
          # UN[k] <- gsub(BLOCK[2],pattern="\\s+",replacement = "") |>
          #   str_split( pattern =  "[(?<=\\(+)][(?=\\)+)]") |>
          #   unlist() |> last() |>
          #   str_match(pattern = ParseOptions$PATTERN_UN ) |>
          #   tolower()

          UN[k] <- gsub(RAW[HL+2],pattern="\\s+",replacement = "") |>
            str_split( pattern =  "[(?<=\\(+)][(?=\\)+)]") |>
            unlist() |> last() |>
            str_match(pattern = ParseOptions$PATTERN_UN ) |>
            tolower() |> as.vector()
        }

        if(TypeID %in% c("TRB","TRC")){
          UN[k] <-  str_match_all(string =HEADER,pattern =  ParseOptions$PATTERN_UN) |>
            unlist()  |>  str_remove_all(pattern = "\\s*") |> tolower()
        }
        if(length(UN[k])!=1){return(NULL)}

        if(TypeID %in% c("TRZ","TRA")){
          # AUX <- read_table(BLOCK,skip=2,col_names =  FALSE ) |> as.data.table()
          LAST_COL <- colnames(BLOCK)[ncol(BLOCK)]
          AUX <- BLOCK[,..LAST_COL]
        }
        if(TypeID %in% c("TRB","TRC")){
          # AUX <- read_table(BLOCK,skip=0,col_names =  FALSE ) |> as.vector()
          # AUX <- read_table(BLOCK,skip=0,col_names =  FALSE ) |> as.data.table()
          AUX <- BLOCK
        }
        TSL <- c(TSL,list(AUX))
        # cat(sprintf("Done\n"))
      }
      if(all(UN==UN[1])){
        UN <- unique(UN)
      } else {return(NULL)}

      if(grepl(UN,pattern = "[///+]")){
        UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
      }
      if(any(dt!=dt[1])){return(NULL)}

      dt <- unique(dt)
      NP <- max(lengths(TSL))



      for(k in seq_len(length.out = nrow(FT))){
        if(length(TSL[[k]])<NP){length(TSL[[k]]) <- NP}
      }
      # TS <- cbind(
      #   seq(from=0,length.out=NP,by=dt),
      #   as.data.table(TSL) |> setnafill(fill = 0))
      # colnames(TS) <- c("T",DIR)
      TS <- as.data.table(TSL) |> setnafill(fill = 0)
      colnames(TS) <- DIR


      # Parse HEADER: OwnerEventSN (TRA only)
      # if(is.na(OwnerEventSN) && TypeID %in% c("TRA","TRZ")){
      # if(TypeID %in% c("TRA","TRZ")){
      #   AUX <- str_match_all(HEADER,pattern = "(?<=Event\\sdate\\:).*") |> unlist() |> str_trim()
      #   MMDDYY <- str_split(AUX,"/",simplify = TRUE)
      #   if(length(MMDDYY)==3){
      #     YY <- MMDDYY[1] |> str_trim()
      #     DD <- MMDDYY[3] |> str_trim()
      #     if(nchar(DD)==1) DD <- paste0(rep("0",2-length(DD)),DD)
      #     MM <- MMDDYY[2] |> str_trim()
      #     if(nchar(MM)==1) MM <- paste0(rep("0",2-length(MM)),MM)}
      #   else {
      #     YY <- "0000"
      #     MM <- "00"
      #     DD <- "00"}
      #   AUX <- str_match_all(HEADER,pattern = "(?<=Event\\sorigin\\stime\\:).*") |> unlist() |> str_trim()
      #   AUX <- str_split(AUX,pattern ="[.]",simplify = TRUE)[1]
      #   HHMMSS <- str_split(AUX,pattern = "[:]") |> unlist() |> str_trim()
      #   if(length(HHMMSS)==3){
      #     hh <- HHMMSS[1]
      #     if(nchar(hh)==1) hh <- paste0(rep("0",2-length(hh)),hh)
      #     mm <- HHMMSS[2]
      #     if(nchar(mm)==1) mm <- paste0(rep("0",2-length(mm)),mm)
      #     ss <- HHMMSS[3]
      #     if(nchar(ss)==1) ss <- paste0(rep("0",2-length(ss)),ss)
      #
      #   }
      #   else {
      #     hh <- "00"
      #     mm <- "00"
      #     ss <- "00"}
      #   OwnerEventSN <- paste0(YY,DD,MM,hh,mm,ss)
      # }


      OK <<-TRUE
    # },
    # error = function(e){OK <<-FALSE},
    # warning = function(e){OK <<-FALSE}
  # )
  if(OK){
    rawGMSR <- list(
      isRawRecord = TRUE,
      isRegularRecord= FALSE,
      RecordSN = RecordSN,
      OwnerRecordSN = OwnerRecordSN,
      OwnerStationSN = OwnerStationSN,
      # OwnerEventSN = OwnerEventSN,
      Header = gsub('[^\x20-\x7E]', '', HEADER),
      Units = UN,
      dt = dt,
      AT = TS,
      VT = NA,
      DT = NA)
  }

  return(rawGMSR)
}

.readAT <- function(FT=data.table()){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  rawGMSR <- NULL
  OK <- is.data.table(FT) && FT$TypeID[1] %in% c("AT2")
  stopifnot(OK)
  # tryCatch(
  #   expr = {
      DIR <- vector(mode = "character",length = nrow(FT))
      dt<- vector(mode = "double",length = nrow(FT))
      NP <- vector(mode = "double",length = nrow(FT))
      UN <- vector(mode = "character",length = nrow(FT))
      TSL <- list()

      for(k in seq_len(length.out = nrow(FT))){
        cat(sprintf(" %s ",FT$OwnerFileSN[k]))
        TypeID <- FT$TypeID[k]
        OwnerRecordSN <- FT$OwnerRecordSN[k]
        OwnerStationSN <- FT$OwnerStationSN[k]
        # OwnerEventSN <- FT$OwnerEventSN[k]
        RecordSN <- FT$RecordSN[k]

        ParseOptions <- .setParseOptions(TypeID)
        stopifnot(!is.null(ParseOptions))
        OwnerFilename <- FT$OwnerFilename[k]

        RAW <- read_lines(OwnerFilename,skip_empty_rows = FALSE)
        RAW <- RAW[RAW!=""]
        # HL <- ParseOptions$HEADER_LINES
        HL <- NULL
        HL <- grep(x=RAW,pattern = "NPTS")
        stopifnot(length(HL)==1 & HL>0)
        HEADER <- RAW[c(1:HL)]

        DIR[k] <- str_split(HEADER[2],pattern = "[,]") |> unlist() |> last() |>
          str_remove_all(pattern = "\\s")
        stopifnot(!is.na(DIR[k]))

        dt[k] <- str_match_all(string =HEADER,pattern =  ParseOptions$PATTERN_DT) |>
          unlist()  |>  str_remove_all(pattern = "\\s*") |> as.numeric()
        stopifnot(length(dt[k])==1 & dt[k]>0)

        UN[k] <- str_match_all(string =HEADER, pattern =  ParseOptions$PATTERN_UN) |>
          unlist()  |> str_remove_all(pattern = "\\s*") |> tolower()
        stopifnot(length(UN[k])==1)

        # Fix 123.45-234.11 -> 123.45 -123.11
        BLOCK <- RAW[-c(1:HL)] |>
          str_replace_all(pattern = "(?<=[[:digit:]])([-])(?=[[:digit:]])",replacement =" -")
        AT <- fread(text = BLOCK, fill = TRUE, data.table = TRUE) |>
          transpose() |> unlist() |> as.vector()
        TSL <- c(TSL,list(AT))

        # AUX <- fread(text = BLOCK, fill = TRUE, data.table = TRUE) |>
        #   unlist() |> as.vector() #as.data.table()
        # TSL <- c(TSL,list(AUX))
        # cat(sprintf("Done\n"))
      }
      stopifnot(all(UN==UN[1]))
      UN <- unique(UN)
      if(grepl(UN,pattern = "[///+]")){
        UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
      }
      stopifnot(all(dt==dt[1]))
      dt <- unique(dt)
      NP <- max(lengths(TSL))

      for(k in seq_len(length.out = nrow(FT))){
        if(length(TSL[[k]])<NP){length(TSL[[k]]) <- NP}
      }
      # TS <- cbind(
      #   seq(from=0,length.out=NP,by=dt),
      #   as.data.table(TSL) |> setnafill(fill = 0))
      # colnames(TS) <- c("T",DIR)
      TS <- as.data.table(TSL) |> setnafill(fill = 0)
      colnames(TS) <- DIR



      # Parse HEADER: OwnerEventSN
      # if(is.na(OwnerEventSN)){
      #   AUX <- str_split(HEADER[2], pattern = "[,]", simplify = TRUE)
      #   MMDDYY <- str_split(AUX[2],"/",simplify = TRUE)
      #   YY <- MMDDYY[3] |> str_trim()
      #   DD <- MMDDYY[2] |> str_trim()
      #   if(nchar(DD)==1) DD <- paste0(rep("0",2-length(DD)),DD)
      #   MM <- MMDDYY[1] |> str_trim()
      #   if(nchar(MM)==1) MM <- paste0(rep("0",2-length(MM)),MM)
      #   OwnerEventSN <- paste0(YY,DD,MM,"000000")
      # }

      OK <<-TRUE
  #   },
  #   error = function(e){OK <<-FALSE},
  #   warning = function(e){OK <<-FALSE}
  # )
  if(OK){
    rawGMSR <- list(
      isRawRecord = TRUE,
      isRegularRecord= FALSE,
      RecordSN = RecordSN,
      OwnerRecordSN = OwnerRecordSN,
      OwnerStationSN = OwnerStationSN,
      # OwnerEventSN = OwnerEventSN,
      Header = gsub('[^\x20-\x7E]', '', HEADER),
      Units = UN,
      dt = dt,
      AT = TS,
      VT = NA,
      DT = NA)}
  return(rawGMSR)
}

.readAC <- function(FT=data.table()){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  rawGMSR <- NULL
  OK <- is.data.table(FT) && nrow(FT)==1 && FT$TypeID[1] %in% c("ACA","ACB","LIS")
  stopifnot(OK)
  tryCatch(
    expr = {
      TypeID <- FT$TypeID
      OwnerRecordSN <- FT$OwnerRecordSN
      OwnerStationSN <- FT$OwnerStationSN
      RecordSN <- FT$RecordSN

      OwnerFilename <- FT$OwnerFilename
      cat(sprintf(" %s ",basename(OwnerFilename)))
      # cat(sprintf("> Parsing %s file: %s ...",TypeID,OwnerFilename))
      ParseOptions <- .setParseOptions(TypeID)
      stopifnot(!is.null(ParseOptions))
      HL <- NULL
      RAW <- read_lines(OwnerFilename,skip_empty_rows = TRUE)
      if(TypeID == "ACB"){
        HL <- RAW |> grep(pattern = "\\s+[T]\\s+[EW]")-1
      }
      if(TypeID == "ACA"){
        HL <- RAW |> grep(pattern = "\\s+[Z]\\s+[N]")-1
      }
      if(TypeID == "LIS"){
        HL <- RAW |> grep(pattern = "\\={2,}DATA")
      }
      stopifnot(length(HL)>0 & HL>0)
      HEADER <- RAW[c(1:HL)] |> gsub(pattern = "\\n",replacement = ".") |> as.vector()
      dt <- str_match_all(string = HEADER,pattern =  ParseOptions$PATTERN_DT) |>
        unlist() |>  str_remove_all(pattern = "\\:") |>  str_remove_all(pattern = "\\s*")
      stopifnot(length(dt)==1 & dt>0)
      dt <- dt |> as.numeric()

      if(TypeID %in% c("ACA","ACB")){ dt <- 1/dt}
      if(TypeID %in% c("ACA","ACB")){
        UN <- str_match_all(string = HEADER, pattern =  ParseOptions$PATTERN_UN) |>
          unlist() |> str_remove_all(pattern = "\\:") |> str_remove_all(pattern = "\\s*")
        stopifnot(length(UN)==1)
        UN <- UN[UN!=""]  |> tolower()
        if(grepl(UN,pattern = "[///+]")){
          UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
        }
      }
      if(TypeID=="LIS"){UN <- "cm"}

      BLOCK <- fread(text=RAW[HL+1:length(RAW)],fill = TRUE) |> na.omit() # RAW[-c(1:HL)]
      #Check if the record has Column Names
      # stopifnot(!.allNumbers(BLOCK[1]))
      TS <- BLOCK # read_table(BLOCK,col_names =  TRUE) |> as.data.table()
      NP <- nrow(TS)
      if(TypeID == "ACB")  {
        # remove time
        TS[,T:=NULL]
        # TS <- cbind(TS,data.table(T=seq(from=0,length.out=NP,by=dt)))
      }
      OK <<-TRUE
    },
    error = function(e){OK <<-FALSE},
    warning = function(e){OK <<-FALSE}
  )
  if(OK){
    rawGMSR <- list(
      isRawRecord = TRUE,
      isRegularRecord= FALSE,
      RecordSN = RecordSN,
      OwnerRecordSN = OwnerRecordSN,
      OwnerStationSN = OwnerStationSN,
      # OwnerEventSN = OwnerEventSN,
      Header = gsub('[^\x20-\x7E]', '', HEADER),
      Units = UN,
      dt = dt,
      AT = TS,
      VT = NA,
      DT = NA)}
  return(rawGMSR)
}

.readV2 <- function(FT=data.table()){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  rawGMSR <- NULL
  OK <- is.data.table(FT) & nrow(FT)==1 & FT$TypeID[1] %in% c("V2A","V2C")
  stopifnot(OK)
  # tryCatch(
  #   expr = {
  #
      TypeID <-FT$TypeID
      OwnerRecordSN <- FT$OwnerRecordSN
      OwnerStationSN <- FT$OwnerStationSN
      # OwnerEventSN <- FT$OwnerEventSN
      RecordSN <- FT$RecordSN
      OwnerFilename <- FT$OwnerFilename
      ParseOptions <- .setParseOptions(TypeID)
      stopifnot(!is.null(ParseOptions))
      cat(sprintf(" %s ",OwnerFilename))
      RAW <- read_lines(OwnerFilename,skip_empty_rows = FALSE)
      RAW <- RAW[RAW!=""]
      HL <- NULL
      if(TypeID=="V2A"){
        HL <- grep(x=RAW,pattern = "(?i)Displacement\\:")|> first()
        HL <- HL+10
      }

      if(TypeID=="V2C"){
        HL <- grep(x=RAW,pattern = "POINTS\\sOF\\sACCEL") |> first()
      }
      stopifnot(length(HL)==1 & HL>0)

      HEADER <- RAW[c(1:HL)] |> as.vector() #as.data.table()
      DIR <- str_match_all(string = RAW,pattern =  ParseOptions$PATTERN_DIR) |>
        unlist() |>  str_remove_all(pattern = "\\s") |> as.vector()
      if(length(DIR)==4){
        DIR <- DIR[-4]}
      # stopifnot(length(DIR)==3)
      if(length(DIR)!=3) return(NULL)

      dt <- str_match_all(string = HEADER,pattern =  ParseOptions$PATTERN_DT) |>
        unlist() |>  str_remove_all(pattern = "\\s") |> unique()
      stopifnot(length(dt)==1)
      dt <- dt[dt!=""]  |> as.numeric()

      NP <- str_match_all(string =RAW,pattern =  ParseOptions$PATTERN_NP) |>
        unlist() |> str_remove_all(pattern = "\\s")
      NP <- NP[NP!=""]
      if(length(NP)==4){
        NP <- NP[-4]}
      # stopifnot(length(NP)==3 )
      if(length(NP)!=3) return(NULL)
      NP <- NP |> as.numeric()
      stopifnot(all(NP>0))

      UN <- str_match_all(string =HEADER, pattern =  ParseOptions$PATTERN_UN) |>
        unlist() |> str_remove_all(pattern = "\\s") |> last() |> tolower()
      stopifnot(length(UN)==1)
      if(grepl(UN,pattern = "[///+]")){
        UN <- (str_split(UN, pattern = "[///+]") |> unlist())[1]
      }
      if(TypeID=="V2C") {
        SHIFT <- 1# c(1,1,1,0)
        IDX <- grep(x=RAW,pattern = ParseOptions$PATTERN_TS_START)+SHIFT
        IDX <- c(IDX,length(RAW)+1)
      }

      if(TypeID=="V2A") {
        SHIFT <- HL#  c(HL,HL,HL,0)
        IDX <- grep(x=RAW,pattern = ParseOptions$PATTERN_BLOCK_START)+SHIFT
        IDX <- c(IDX,length(RAW)+1)
      }
      stopifnot(length(IDX)==4)
      TSL <- list()
      NCOLS <- ParseOptions$NCOLS
      for(j in seq_along(DIR)){
        BLOCK <- RAW[IDX[j]:(IDX[j+1]-1)]
        NR <- floor(NP[j]/NCOLS)
        # Fix 123.45-234.11 -> 123.45 -123.11
        BLOCK <- BLOCK[1:NR] |>
          str_replace_all(pattern = "(?<=[[:digit:]])([-])(?=[[:digit:]])",replacement =" -")
        AT <- fread(text = BLOCK, fill = TRUE, data.table = TRUE) |>
          transpose() |> unlist() |> as.vector()
        TSL <- c(TSL,list(AT))
      }

      NP <- max(lengths(TSL))+NCOLS # Aqui NP ya es un escalar
      for(k in seq_along(DIR)){
        if(length(TSL[[k]])<NP){length(TSL[[k]]) <- NP}
      }
      # TS <- cbind(
      #   seq(from=0,length.out=NP,by=dt),
      #   as.data.table(TSL) |> setnafill(fill = 0))
      # colnames(TS) <- c("T",DIR)

      TS <- as.data.table(TSL) |> setnafill(fill = 0)
      colnames(TS) <- DIR
      # Parse HEADER: OwnerEventSN
      # if(is.na(OwnerEventSN) && TypeID=="V2A"){
      #   AUX <- str_split(HEADER[8],pattern = "[\\s+]",simplify = TRUE)
      #   YYMMDD <- AUX[AUX!=""]
      #   if(length(YYMMDD)>=4){
      #     YY <- YYMMDD[1] |> str_trim()
      #     DD <- YYMMDD[3] |> str_trim()
      #     if(nchar(DD)==1){DD <- paste0(rep("0",2-length(DD)),DD)}
      #     MM <- match(YYMMDD[2] |> str_trim(),month.name) |> as.character()
      #     if(nchar(MM)==1) {MM <- paste0(rep("0",2-length(MM)),MM)}
      #   }
      #   else {
      #     YY <- "0000"
      #     MM <- "00"
      #     DD <- "00"}
      #   HHMMSS <- paste0(AUX[4],"00")
      #   OwnerEventSN <- paste0(YY,MM,DD,HHMMSS)
      # }

      # if(is.na(OwnerEventSN) && TypeID=="V2C"){
      #   AUX <- str_split(HEADER[4],pattern = "[(GMT)]",simplify = TRUE)[1]
      #   YYMMDD <- str_split(AUX,pattern = "[\\s+]",simplify = TRUE)[1] |> str_split(pattern = "[/]",simplify = TRUE)
      #
      #   if(length(YYMMDD)==3){
      #     YY <- YYMMDD[3] |> str_trim()
      #     DD <- YYMMDD[1] |> str_trim()
      #     if(nchar(DD)==1){DD <- paste0(rep("0",2-length(DD)),DD)}
      #     MM <- YYMMDD[2] |> str_trim()
      #     if(nchar(MM)==1) {MM <- paste0(rep("0",2-length(MM)),MM)}
      #   }
      #   else {
      #     YY <- "0000"
      #     MM <- "00"
      #     DD <- "00"}
      #
      #   HHMMSS <- str_split(AUX,pattern = "[\\s+]",simplify = TRUE)[2] |> str_split(pattern = "[:]",simplify = TRUE)
      #   if(length(HHMMSS)==3){
      #     hh <- HHMMSS[1]
      #     if(nchar(hh)==1) hh <- paste0(rep("0",2-length(hh)),hh)
      #     mm <- HHMMSS[2]
      #     if(nchar(mm)==1) mm <- paste0(rep("0",2-length(mm)),mm)
      #     ss <- HHMMSS[3]
      #     if(nchar(ss)==1) ss <- paste0(rep("0",2-length(ss)),ss)
      #
      #   }
      #   else {
      #     hh <- "00"
      #     mm <- "00"
      #     ss <- "00"}
      #   OwnerEventSN <- paste0(YY,DD,MM,hh,mm,ss)
      #
      # }
      OK <<-TRUE
  #   },
  #   error = function(e){OK <<-FALSE},
  #   warning = function(e){OK <<-FALSE}
  # )
  if(OK){
    rawGMSR <- list(
      isRawRecord = TRUE,
      isRegularRecord= FALSE,
      RecordSN = RecordSN,
      OwnerRecordSN = OwnerRecordSN,
      OwnerStationSN = OwnerStationSN,
      # OwnerEventSN = OwnerEventSN,
      Header = gsub('[^\x20-\x7E]', '', HEADER),
      Units = UN,
      dt = dt,
      AT = TS,
      VT = NA,
      DT = NA)}
  return(rawGMSR)
}

.parseOwnerHeader<- function(TypeID, HEADER=NULL){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  stopifnot(toupper(TypeID) %in% c("TRZ","TRA","TRB","TRC","V2A","ACA","ACB","LIS","V2C","AT2") && !is.null(HEADER))
  HEADER <- gsub('[^\x20-\x7E]', '', HEADER)
  # Tag V2A data ------------------------------------------------------------
  if(TypeID =="V2A" ){
    StationName <- HEADER[3]
    TOKEN <- str_split(string = HEADER[2],pattern = "\\s+") |> unlist()
    StationSN <- TOKEN[2]
    StationLatitude <- paste0(TOKEN[3],".",TOKEN[4],".",TOKEN[5])|> trimws()
    StationLongitude <- paste0(TOKEN[6],".",TOKEN[7],".",TOKEN[8])|> trimws()
    EventName <- HEADER[7]
    TOKEN <- str_split(string = HEADER[8],pattern = "\\s+") |> unlist()
    EventDate <- paste0(TOKEN[1],".",TOKEN[2],".",TOKEN[3])|> trimws()
    EventTime <- paste0(TOKEN[4],".",TOKEN[5])|> trimws()
    TOKEN <- str_split(string = HEADER[9],pattern = "\\s+") |> unlist()|> trimws()
    EpicenterLatitude <- paste0(TOKEN[2],".",TOKEN[3],".",TOKEN[4])|> trimws()
    EpicenterLongitude <- paste0(TOKEN[5],".",TOKEN[6],".",TOKEN[7])|> trimws()
    EpicenterDistance <- TOKEN[11]|> trimws()
    HypocenterDepth <- TOKEN[13]|> trimws()
    Magnitude <- TOKEN[16]|> trimws()
    MagnitudeType <- TOKEN[15]|> trimws()
  }

  # Tag V2C data ------------------------------------------------------------
  if(TypeID =="V2C" ){
    StationName <- NA
    StationSN <- NA
    StationLatitude <- NA
    StationLongitude <- NA
    EventName <- NA
    EventDate <- NA
    EventTime <- NA
    EpicenterLatitude <- NA
    EpicenterLongitude <- NA
    EpicenterDistance <- NA
    HypocenterDepth <- NA
    Magnitude <- NA
    MagnitudeType <- NA
  }
  # Tag ACA data ------------------------------------------------------------

  if(TypeID =="ACA" ){
    POS <- grep(HEADER,pattern = "1. ES")
    StationName <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationSN <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLatitude <- str_split(HEADER[POS+3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLongitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    POS <- grep(HEADER,pattern = "2. SISMO")
    EventName <- NA
    EventDate <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EventTime <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLatitude <- str_split(HEADER[POS+3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLongitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterDistance <- str_split(HEADER[POS+7],pattern="[:]",simplify = TRUE)[2]|> trimws()
    HypocenterDepth <- str_split(HEADER[POS+5],pattern="[:]",simplify = TRUE)[2]|> trimws()
    Magnitude <- str_extract(HEADER[POS+6],pattern = "(?<=\\:).*[0-9]+\\.[0-9]+") #|> unlist() |> paste(collapse="")
    MagnitudeType <- NA

  }
  # Tag ACB data ------------------------------------------------------------

  if(TypeID =="ACB" ){
    POS <- grep(HEADER,pattern = "THE SEISMIC STATION")
    StationName <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationSN <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLatitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLongitude <- str_split(HEADER[POS+5],pattern="[:]",simplify = TRUE)[2]|> trimws()
    POS <- grep(HEADER,pattern = "THE EARTHQUAKE")
    EventName <- NA
    EventDate <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    TOKEN <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)|> trimws()
    EventTime <- paste0(TOKEN[2],":",TOKEN[3],":",TOKEN[4])
    EpicenterLatitude <- str_split(HEADER[POS+3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLongitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterDistance <- NA
    HypocenterDepth <- str_split(HEADER[POS+5],pattern="[:]",simplify = TRUE)[2]|> trimws()
    TOKEN <- str_split(HEADER[POS+6],pattern="[:]",simplify = TRUE)[2]|> trimws()
    Magnitude <- str_split(TOKEN,pattern="\\s+",simplify = TRUE)[1]|> trimws()
    MagnitudeType <- str_split(TOKEN,pattern="\\s+",simplify = TRUE)[2]|> trimws()

  }
  # Tag LIS data ------------------------------------------------------------
  if(TypeID =="LIS" ){
    StationName <-  str_extract(HEADER[4], pattern = "(?<=Station\\s(?i)name).*")  |> trimws()
    StationSN <- str_split(HEADER[13],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLatitude <- str_split(HEADER[16],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLongitude <- str_split(HEADER[17],pattern="[:]",simplify = TRUE)[2]|> trimws()

    TOKEN <- (str_split(HEADER[7],pattern="Event\\sDate\\:",simplify = TRUE) |> unlist())[2]|> trimws()
    EventName <- str_extract(HEADER[3], pattern = "(?<=^Epicenter\\s).*$")  |> trimws()
    EventDate <- str_split(TOKEN,pattern = " ",simplify = TRUE)[1]
    EventTime <- str_split(TOKEN,pattern = " ",simplify = TRUE)[2]
    EpicenterLatitude <- str_split(HEADER[8],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLongitude <- str_split(HEADER[9],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterDistance <- str_split(HEADER[21],pattern="[:]",simplify = TRUE)[2]|> trimws()
    HypocenterDepth <- str_split(HEADER[22],pattern="[:]",simplify = TRUE)[2]|> trimws()
    Magnitude <- str_split(HEADER[11],pattern="[:]",simplify = TRUE)[2]|> trimws()
    MagnitudeType <- "Mw"

  }
  # Tag TRB data ------------------------------------------------------------
  if(TypeID =="TRB" ){
    StationName <- NA
    StationSN <- str_extract(HEADER[4],pattern = "(?<=(?i)Estacion\\:).*(?=Componente\\:)") |> trimws()
    StationLatitude <- str_extract(HEADER[5],pattern = "(?<=(?i)Latitud\\:).*(?=Longitud\\:)") |> trimws()
    StationLongitude <- str_extract(HEADER[5],pattern = "(?<=(?i)Longitud\\:).*") |> trimws()
    TOKEN <- str_extract(HEADER[1],pattern="(?<=Origen\\:).*")|> trimws()
    EventName <- NA
    EventDate <- str_split(TOKEN,pattern="T",simplify = TRUE)[1]|> trimws()
    EventTime <- str_split(TOKEN,pattern="T",simplify = TRUE)[2]|> trimws()
    EpicenterLatitude <- NA
    EpicenterLongitude <- NA
    EpicenterDistance <- NA
    HypocenterDepth <- NA
    Magnitude <- NA
    MagnitudeType <- NA
  }

  # Tag TRC data ------------------------------------------------------------
  if(TypeID =="TRC" ){
    StationName <- str_split(HEADER[16],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationSN <- str_split(HEADER[15],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLatitude <- str_split(HEADER[17],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLongitude <-  str_split(HEADER[18],pattern="[:]",simplify = TRUE)[2]|> trimws()

    EventName <- str_split(HEADER[1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EventDate <- str_split(HEADER[3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EventTime <- str_split(HEADER[4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLatitude <- str_split(HEADER[5],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLongitude <- str_split(HEADER[6],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterDistance <- str_split(HEADER[25],pattern="[:]",simplify = TRUE)[2]|> trimws()
    HypocenterDepth <- str_split(HEADER[7],pattern="[:]",simplify = TRUE)[2]|> trimws()
    Magnitude <- str_split(HEADER[9],pattern="[:]",simplify = TRUE)[2]|> trimws()
    MagnitudeType <- "Mw"


  }
  # Tag TRA data ------------------------------------------------------------
  if(TypeID =="TRA" ){
    POS <- grep(HEADER,pattern = "STATION INFORMATION")
    StationName <- NA
    StationSN <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLatitude <- str_split(HEADER[POS+3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLongitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    POS <- grep(HEADER,pattern = "EVENT INFORMATION")
    EventName <- NA
    EventDate <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EventTime <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLatitude <- str_split(HEADER[POS+3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLongitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterDistance <- str_split(HEADER[POS+7],pattern="[:]",simplify = TRUE)[2]|> trimws()
    HypocenterDepth <- str_split(HEADER[POS+5],pattern="[:]",simplify = TRUE)[2]|> trimws()#str_match(HEADER[POS+5],pattern = "(?<=\\:).*") |> unlist() |> paste(collapse="")
    Magnitude <- str_split(HEADER[POS+6],pattern="[:]",simplify = TRUE)[2]|> trimws()
    MagnitudeType <- NA

  }
  # Tag TRZ data ------------------------------------------------------------
  if(TypeID =="TRZ" ){
    POS <- grep(HEADER,pattern = "STATION INFORMATION")
    StationSN <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationName <- NA
    StationLatitude <- str_split(HEADER[POS+5],pattern="[:]",simplify = TRUE)[2]|> trimws()
    StationLongitude <- str_split(HEADER[POS+6],pattern="[:]",simplify = TRUE)[2]|> trimws()
    POS <- grep(HEADER,pattern = "EVENT INFORMATION")
    EventName <- NA
    EventDate <- str_split(HEADER[POS+1],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EventTime <- str_split(HEADER[POS+2],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLatitude <- str_split(HEADER[POS+3],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterLongitude <- str_split(HEADER[POS+4],pattern="[:]",simplify = TRUE)[2]|> trimws()
    EpicenterDistance <- str_split(HEADER[POS+7],pattern="[:]",simplify = TRUE)[2]|> trimws()
    HypocenterDepth <- str_split(HEADER[POS+5],pattern="[:]",simplify = TRUE)[2]|> trimws()
    Magnitude <- str_split(HEADER[POS+6],pattern="[:]",simplify = TRUE)[2]|> trimws()
    MagnitudeType <- NA

  }
  # Tag AT2 data ------------------------------------------------------------
  # Parkfield-02 CA, 9/28/2004, PARKFIELD - UPSAR 08, 360
  if(TypeID =="AT2" ){
    TOKEN <- str_split(HEADER[2],pattern=",",simplify = TRUE)
    StationSN <- NA
    StationName <-  TOKEN[3]
    StationLatitude <- NA
    StationLongitude <- NA
    EventName <- TOKEN[1]
    EventDate <- TOKEN[2]
    EventTime <- NA
    EpicenterLatitude <- NA
    EpicenterLongitude <- NA
    EpicenterDistance <- NA
    HypocenterDepth <- NA
    Magnitude <- NA
    MagnitudeType <- NA
  }
  # return  ------------------------------------------------------------
  return(data.table(
    OwnerStationSN = StationSN,
    StationName =  StationName,
    StationLatitude = StationLatitude,
    StationLongitude = StationLongitude,
    EventName = EventName,
    EventDate = EventDate,
    EventTime = EventTime,
    EpicenterLatitude = EpicenterLatitude,
    EpicenterLongitude = EpicenterLongitude,
    HypocenterDepth = HypocenterDepth,
    Magnitude = Magnitude,
    MagnitudeType = MagnitudeType,
    EpicenterDistance = EpicenterDistance
  ))
}

.buildOwnerHeaderTables <- function(TypeID,GMSP=NULL,OVERRIDE=FALSE){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  #Local Copy
  TID <- TypeID
  PATH <- file.path(SourceFolder,"Rds",TypeID)
  # Get EXT
  EXT <- "*.Rds"
  SET <- data.table(
    RecordFilename=list.files(
      path = PATH,
      pattern = EXT,
      recursive = TRUE,
      ignore.case=TRUE,
      include.dirs = FALSE)
  )
  NR <- nrow(SET)
  cat(sprintf("> Processing %d Owner Headers in %s database ...",NR,TypeID))
  # RDSfile_HT <- file.path(IndexFolder,"Rds",paste0(TID,".HeaderTable.Rds"))
  RDSfile_HT <- file.path(IndexFolder,paste0(TypeID,".HeaderTable.Rds"))
  if(!file.exists(RDSfile_HT) || OVERRIDE==TRUE){
    HT <- NULL
    for(i in seq_len(length.out = NR)){
      AUX <- data.table()
      RDSfile_GMSR <- file.path(SourceFolder,"Rds",TypeID,SET$RecordFilename[i])
      GMSR <- readRDS(RDSfile_GMSR)
      AUX <- .parseOwnerHeader(TypeID=TID,HEADER=GMSR$Header)
      AUX[,TypeID:=TID]
      AUX[,OwnerRecordSN:=GMSR$OwnerRecordSN]
      AUX[,RecordSN:=GMSR$RecordSN]
      # AUX[,OwnerEventSN:=GMSR$raw$OwnerEventSN]
      # AUX[,OwnerStationSN:=GMSR$raw$OwnerStationSN]
      HT <- rbindlist(list(HT,AUX),use.names = TRUE)
    }


    saveRDS(compress=FALSE, object=HT, file = RDSfile_HT)
    cat(sprintf("> Done!.\n"))
  }
  else {
    cat(sprintf("* Skip (Header Table already exists)\n"))
    HT <- readRDS(RDSfile_HT)
  }
  return(HT)
}

.allNumbers <- function(TEXT){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  AUX <- TEXT |> str_split(pattern = "\\s") |> unlist()
  return(AUX[AUX!=""] |> grepl(pattern = "[[:digit:]]") |> all())
}

## -- 2. Ground-Motion Intensity Measures(GMIM) -------------------------------------------------------
# Copyright SRK (C) 2010-2021
# Author: A.Verri K. (averri@srk.com.ar)
# Maintainer: A.Verri K.
# Contributors: M. Balbi, P. Barbieri, J. Mussat
# Date: 26-05-2021
#

.buildOwnerIntensityTable <- function(TypeID,GMSP=NULL,MULTI=FALSE, OVERRIDE=FALSE){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  stopifnot(!is.null(GMSP) && toupper(TypeID) %in% c("TRZ","TRA","TRB","TRC","V2A","ACA","ACB","LIS","V2C","AT2"))
  # Build FileTable - ------------------------------------------------------
  # Set Input/Output Paths
  PATH <- file.path(SourceFolder,"Rds",TypeID)
  cat(sprintf("> Building Index #%s ... from %s",TypeID,PATH))
  EXT <- "*.Rds"
  SET <- data.table(
    SourceFilename=list.files(
      path = PATH,
      pattern = EXT,
      recursive = TRUE,
      ignore.case=TRUE,
      include.dirs = FALSE)
  )
  NR <- nrow(SET)
  if(MULTI==TRUE){SET <- SET[sample(1:NR)]}
  # Build Intensity - ------------------------------------------------------
  # RDSfile_IM <- file.path(IndexFolder,"Rds",paste0(TypeID,".IntensityTable.Rds"))
  RDSfile_IM <- file.path(IndexFolder,paste0(TypeID,".OwnerIntensityTable.Rds"))
  if(!file.exists(RDSfile_IM) || OVERRIDE==TRUE){
    IM <- NULL
    NR <- nrow(SET)
    for(i in seq_len(length.out = NR)){
      # Get Record --------------------------------
      cat(sprintf("> Processing record #%s (%d/%d) ...\n",SET$SourceFilename[i],i,NR))
      OK <- TRUE
      RDSfile_GMSR <- file.path(PATH,SET$SourceFilename[i])
      cat(sprintf("     Reading source record (%s)...",basename(RDSfile_GMSR))) # File always exists
      GMSR <- readRDS(RDSfile_GMSR)
      cat(sprintf(" Done!.\n"))
      # Get Intensity Table --------------------------------
      cat(sprintf("     Getting Intensity Measures (IM)..."))

      ##  get raw Record Intensity
      tryCatch(
        expr = {
          AUX <- .getIntensity(a=GMSR$AT,dt=GMSR$dt,UN=GMSR$Units,GMSP=GMSP)
          OK <- !is.null(AUX)
        },
        error = function(e){
          cat(sprintf("* ERROR: %s *",e))
          OK <<- FALSE},
        warning = function(e){
          cat(sprintf("* WARNING: %s *",e))
          OK <<- FALSE}
      )
      if(OK){
        AUX[,TSID:="ATH"]
        AUX[,SiteSN:="I"]
        # AUX[, BlockSN:= NA]
        AUX[,OwnerRecordSN:=GMSR$OwnerRecordSN]
        AUX[,RecordSN:=GMSR$RecordSN]
        AUX[,isRawRecord := TRUE]
        AUX[,isRegularRecord := FALSE]
        AUX[,TypeID:=TypeID]
        IM <- rbindlist(list(IM,AUX),use.names = TRUE)
        cat(sprintf(" Done!.\n"))
      }
      else {cat(sprintf("* Skip (Source Record %s has wrong Units, PGA or NP) \n",SET$RecordSN[i]))}
    }
    cat(sprintf("> Writing Intensity Table #%s ...",RDSfile_IM))
    saveRDS(compress=FALSE, object=IM, file = RDSfile_IM)
    cat(sprintf(" Done!.\n"))
  }
  else {
    cat(sprintf("* Skip (Source Intensity Table already exists)\n"))
    IM <- readRDS(RDSfile_IM)
  }
  return(IM)
}


## -- 3. Ground-Motion Signal Processing (GMSP) ------------------------------------------------------------
# FFT-based integral with Short-Time Fourier Transform
# Copyright SRK (C) 2010-2021
# Author: A. Verri K. (averri@srk.com.ar)
# Contributors: J. Mussat (jmussat@srk.com.ar)
# Maintainer: A.Verri K.
# Date: 26-05-2021

.regularizeOwnerRecords <- function(TypeID,MULTI=TRUE, GMSP=NULL,Amin = 0,Amax=Inf){
  stopifnot(!is.null(GMSP) & !is.null(TypeID))
  TID <- TypeID
  RSN_PROCESSED <- NULL
  Fs <- GMSP$DownFs
  # stopifnot Fs esta en FilterTable
  IM_FILE <- file.path( IndexFolder,paste0(TID,".RegularizedIntensityTable.",Fs,".csv"))
  # Load HeaderIntensityTable - --------------------------------------------
  PATH <- file.path(IndexFolder)
  HT_FILE <- file.path(PATH,"HeaderTable.Rds")
  stopifnot(file.exists(HT_FILE))
  HT <- readRDS(HT_FILE)
  HT <- HT[TypeID==TID] #[!is.na(station.VS30)]
  if(nrow(HT)==0){
    cat(sprintf("> No records available. Exit.\n"))
    return()
  }
  RSN <-  HT$RecordSN  |> unique()
  # Remove records already processed
  if(file.exists(IM_FILE)){
    cat(sprintf("> Reading previous Intensity Table %s...",basename(IM_FILE)))
    AUX <- fread(IM_FILE,yaml = TRUE)
    RSN_PROCESSED <-  AUX$RecordSN  |> unique()
    RSN <- RSN[!(RSN %in% RSN_PROCESSED)]
  }


  # Build FileTable -------------------------------------------------------
  FT <- .buildFileTable(TypeID=TID, MULTI=MULTI)
  SET <- FT[RecordSN %chin% RSN] |> unique(by=c("TypeID","SourceFilename","RecordSN"))
  cat(sprintf("> Found %d regularized records in %s database.\n",nrow(SET),TID))
  NR <- nrow(SET)
  if(NR==0){
    cat(sprintf("> Nothing to do. Exit.\n"))
    return()
  }
  cat(sprintf("> Regularizing %d records ...\n",NR))
  if(MULTI==TRUE){SET <- SET[sample(1:.N)]} # Shuffle records
  # Regularize Record  -------------------------------------------------------
  PATH <- file.path(SourceFolder,"Rds",TID)
  IM <- NULL
  fs <-GMSP$fs
  # for(i in seq_len(length.out = NR)){
  while(nrow(SET)>0){
    i <- 1
    cat(sprintf("> Processing record %d of %d (#%s) from database %s\n",length(RSN_PROCESSED),length(RSN),SET$RecordSN[i],TID))
    # RDSfile_REG <- SET$RecordFilename[i]

    # Get SRC Record  ----
    RDSfile_SRC <- SET$SourceFilename[i]
    GMSR <- list(raw=NULL,reg=NULL)

    # Build TFT  -------------------------------

    VS30 <- HT[RecordSN==RSN[i]]$station.VS30 |> as.numeric()


    # Regularize station Record  --------------------------------
    OUT <- NULL
    tryCatch(
      expr = {
        cat(sprintf("     Reading source file (%s)...",basename(RDSfile_SRC)))
        SRC <- readRDS(RDSfile_SRC)
        cat(sprintf("Done.\n"))
        cat(sprintf("     Normalizing record (%s) ...",SET$RecordSN[i]))
        OUT <- .regularizeRecord(a=SRC$AT,dt=SRC$dt,UN=SRC$Units,GMSP=GMSP,Amin=Amin,Amax=Amax)
        OK <- !is.null(OUT)
        cat(sprintf("Done.\n"))
      },
      error = function(e){
        cat(sprintf("* ERROR: %s *",e))
        OK <<- FALSE},
      warning = function(e){
        cat(sprintf("* WARNING: %s *",e))
        OK <<- FALSE}
    )

    # Get Intensity Table --------------------------------
    IM_SSN <- NULL
    if(OK){
      cat(sprintf("     Getting Intensities ..."))
      TS <- OUT$TS$I #[[j]]
      # Remove unused columns ----
      TS[,W:=NULL]
      COLS_AT <- colnames(TS) |> grep(pattern = paste0("AT","\\."),value = TRUE)
      COLS_VT <- colnames(TS) |> grep(pattern = paste0("VT","\\."),value = TRUE)
      COLS_DT <- colnames(TS) |> grep(pattern = paste0("DT","\\."),value = TRUE)
      AT <- TS[,..COLS_AT]
      VT <- TS[,..COLS_VT]
      DT <- TS[,..COLS_DT]
      # Clean OCID ----
      names(AT) <- str_split(names(AT), pattern = "[[.]]",simplify = TRUE)[,2]
      names(VT) <- str_split(names(DT), pattern = "[[.]]",simplify = TRUE)[,2]
      names(DT) <- str_split(names(DT), pattern = "[[.]]",simplify = TRUE)[,2]
      # Get Intensity ----
      tryCatch(
        expr = {
          AUX <- .getIntensity(a=AT,v=VT,d=DT,dt=OUT$dt,UN=OUT$UN,GMSP=GMSP)
          AUX[,TSID:="ATH"]
          IM_SSN <- rbindlist(list(IM_SSN,AUX))
          IM_SSN[,RecordSN:=SET$RecordSN[i]]
          IM_SSN[,isRawRecord := FALSE]
          IM_SSN[,isRegularRecord := TRUE]
          IM_SSN[,TypeID := TID]
          IM_SSN[,TopLayerID := "station"]
          OK <- !is.null(IM_SSN)
          cat(sprintf("Done.\n"))
        },
        error = function(e){
          cat(sprintf("* ERROR: %s *",e))
          OK <<- FALSE},
        warning = function(e){
          cat(sprintf("* WARNING: %s *",e))
          OK <<- FALSE}
      )
    }
    if(OK){
      # Write CSV ----
      tryCatch(
        expr = {
          cat(sprintf("     Writing Intensity Table ..."))
          fwrite(IM_SSN,file=IM_FILE,yaml = TRUE,append = TRUE)
          # AUX <- fread(IM_FILE,yaml = TRUE)
          RSN_PROCESSED <-  IM_SSN$RecordSN |> unique() #AUX$RecordSN  |> unique()
          RSN <- RSN[!(RSN %in% RSN_PROCESSED)]
          SET <- SET[!(RecordSN %in% RSN_PROCESSED)]
          OK <- TRUE
          cat(sprintf(" Done!.\n"))
        },
        error = function(e){
          cat(sprintf("* ERROR: %s *",e))
          OK <<- FALSE},
        warning = function(e){
          cat(sprintf("* WARNING: %s *",e))
          OK <<- FALSE}
      )
    }
    if(!OK){
      cat(sprintf("* Wrong record. Could not regularize record %s\n",SET$RecordSN[i]))
      SET <- SET[RecordSN!=SET$RecordSN[i]]

    }
    if(nrow(SET)==0){
      cat(sprintf("     Nothing to do. Exit.\n"))
      return()
    }

  }
}

.convoluteOwnerRecords <- function(TypeID,MULTI=FALSE, GMSP=NULL,Amin = 0,Amax=Inf){
  stopifnot(toupper(TypeID) %in% c("TRZ","TRA","TRB","TRC","V2A","ACA","ACB","LIS","V2C","AT2") && !is.null(GMSP))
  if(file.exists(file.path(IndexFolder,"idx.Rds"))){
    idx <- readRDS(file=file.path(IndexFolder,"idx.Rds"))
  } else {
    idx <- NULL
  }
  TID <- TypeID
  # Load HeaderIntensityTable - --------------------------------------------
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"HeaderTable.Rds")
  stopifnot(file.exists(RDSfile))
  HT <- readRDS(RDSfile)
  # Keep Only Records with valid VS30 available ------
  HT <- HT[TypeID==TID][!is.na(station.VS30)]

  if(nrow(HT)==0){
    cat(sprintf("> No records available. Exit.\n"))
    return()
  }
  # Load SiteTable  --------------------------------

  SiteTable <- readRDS(file=file.path(SitesFolder,"SiteTable.Rds"))

  # Load SourceIntensityTable - --------------------------------------------

  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,paste0(TID,".OwnerIntensityTable.Rds"))
  stopifnot(file.exists(RDSfile))
  IM <- readRDS(RDSfile)
  IM <- IM[TypeID==TID]


  RecordSN_VALID <-  HT$RecordSN  |> unique()

  # Keep Only Records with PGA > Amin & PGA < Amax ------
  if(is.null(Amin) & is.null(Amax)) {
    RSN <- IM[RecordSN %chin% RecordSN_VALID,RecordSN] |> unique()
  }
  if(!is.null(Amin) & is.null(Amax)) {
    RSN <- IM[PGA>=Amin & RecordSN %chin% RecordSN_VALID,RecordSN] |> unique()}
  if(is.null(Amin) & !is.null(Amax)) {
    RSN <- IM[PGA<=Amax & RecordSN %chin% RecordSN_VALID,RecordSN] |> unique()
  }
  if(!is.null(Amin) & !is.null(Amax)) {
    RSN <- IM[PGA>=Amin & PGA<=Amax & RecordSN %chin% RecordSN_VALID,RecordSN] |> unique()
  }


  NR_VALID <- length(RSN)
  cat(sprintf("> Found %d Records with valid amplitudes.\n",NR_VALID))
  if(NR_VALID==0){
    cat(sprintf("> Amplitudes out of range. No records available. Exit.\n"))
    return()
  }

  d <- GMSP$D
  # FILE <- file.path(SitesFolder,paste("S2O",Fs,NW,"D",100*d,"Rds",sep="."))
  S2O <-  GMSP$S2O #readRDS(file=FILE)
  # FILE <- file.path(SitesFolder,paste("S2B",Fs,NW,"D",100*d,"Rds",sep="."))
  S2B <-  GMSP$S2B #readRDS(file=FILE)




  # Build FileTable -------------------------------------------------------
  FT <- .buildFileTable(TypeID=TID,MULTI=MULTI)
  SET <- FT[RecordSN %chin% RSN] |> unique(by=c("TypeID","SourceFilename","RecordSN"))
  # cat(sprintf("> Found %d regularized records in %s database.\n",SET[Regularized==TRUE,.N],TID))
  cat(sprintf("> Found %d regularized records in %s database.\n",nrow(SET),TID))
  # if(OVERRIDE==FALSE){SET <- SET[Regularized == FALSE]}
  # SET <- SET[Regularized == FALSE]
  NR <- nrow(SET)
  if(NR==0){
    cat(sprintf("> Nothing to do. Exit.\n"))
    return()
  }
  cat(sprintf("> Regularizing %d records ...\n",NR))
  if(MULTI==TRUE){SET <- SET[sample(1:.N)]} # Shuffle records

  # Build Folder Structure  -------------------------------------------------------
  # BlockSN <-  str_sub(SET$RecordSN,start = 1L, end = 2L)
  # PATH <- file.path(RecordsFolder,"Rds")
  # AUX <- file.path(PATH,BlockSN |> unique())
  # map(AUX,function(x) {dir.create(x,showWarnings = FALSE)})
  # Regularize Record  -------------------------------------------------------

  PATH <- file.path(SourceFolder,"Rds",TID)
  IM <- NULL
  for(i in seq_len(length.out = NR)){
    cat(sprintf("> Processing record #%s from database %s (%d/%d):\n",SET$RecordSN[i],TID,i,NR))
    RDSfile_REG <- SET$RecordFilename[i]

    ## Get SRC Record  ----
    RDSfile_SRC <- SET$SourceFilename[i]
    GMSR <- list(raw=NULL,reg=NULL)
    cat(sprintf("     Reading source record (%s)...",basename(RDSfile_SRC)))
    # GMSR$raw <- readRDS(RDSfile_SRC)
    OUT <- readRDS(RDSfile_SRC)
    cat(sprintf("Done!.\n"))

    ## Build TFT  -------------------------------



    VS30 <- HT[RecordSN==RSN[i]]$station.VS30 |> as.numeric()
    VS30.min <- 0.90*VS30
    VS30.max <- 1.10*VS30
    SID_Valid <- SiteTable[VS30>=VS30.min & VS30<=VS30.max]$SID
    cat(sprintf("     Select %d sites for Vs30 %4.0f m/s...\n",length(SID_Valid),VS30))


    if(!is.null(idx)){
      SID_Processed <- idx[RecordSN==SET$RecordSN[i]]$SiteSN
      cat(sprintf("     Remove %d sites already processed for record %s\n",length(SID_Processed),SET$RecordSN[i]))
      SID_Valid <- SID_Valid[!(SID_Valid %in% SID_Processed)]

    }
    SID=sample(SID_Valid,size=min(length(SID_Valid),GMSP$NSMAX))
    if(length(SID)>0){

      fs <-GMSP$fs


      ## Regularize S2O Record  --------------------------------

      tryCatch(
        expr = {
          cat(sprintf("     Convolute (S2O) record %s through %d sites...",SET$RecordSN[i],length(SID)))
          COLS <- colnames(GMSP$S2O)
          s2o <- GMSP$S2O[, unique(COLS[COLS %in% SID])]  |> as.data.table()#Surface to Outcrop
          s2o[,I:=1]
          OUT <- .regularizeRecord(a=OUT$AT,dt=OUT$dt,UN=OUT$Units,GMSP=GMSP,Amin=Amin,Amax=Amax,TFT=s2o)
          OK <- !is.null(OUT)
          IM_SSN <- NULL
          if(OK){
            cat(sprintf("Done.\n"))
            # GMSR$s2o$isRawRecord = FALSE
            # GMSR$s2o$isRegularRecord = TRUE
            # GMSR$s2o$Units = OUT$UN
            # GMSR$s2o$dt = OUT$dt
            # GMSR$s2o$TS = OUT$TS
            # Get Intensity Table --------------------------------
            cat(sprintf("     Building S2O intensity table ..."))

            # SSN <- names(GMSR$s2o$TS)
            SSN <- names(OUT$TS)
            # ************** debug
            # if("x" %in% SSN){
            #   A <- 1
            # }
            # ************** debug

            for(j in seq_along(SSN)){
              # TS <- GMSR$s2o$TS[[j]]
              TS <- OUT$TS[[j]]
              # Remove unused columns ----
              TS[,W:=NULL]
              COLS_AT <- colnames(TS) |> grep(pattern = paste0("AT","\\."),value = TRUE)
              COLS_VT <- colnames(TS) |> grep(pattern = paste0("VT","\\."),value = TRUE)
              COLS_DT <- colnames(TS) |> grep(pattern = paste0("DT","\\."),value = TRUE)
              AT <- TS[,..COLS_AT]
              VT <- TS[,..COLS_VT]
              DT <- TS[,..COLS_DT]
              # Clean OCID ----
              names(AT) <- str_split(names(AT), pattern = "[[.]]",simplify = TRUE)[,2]
              names(VT) <- str_split(names(DT), pattern = "[[.]]",simplify = TRUE)[,2]
              names(DT) <- str_split(names(DT), pattern = "[[.]]",simplify = TRUE)[,2]
              # Get Intensity ----
              # AUX <- .getIntensity(a=AT,v=VT,d=DT,dt=GMSR$s2o$dt,UN=GMSR$s2o$Units,GMSP=GMSP)
              AUX <- .getIntensity(a=AT,v=VT,d=DT,dt=OUT$dt,UN=OUT$UN,GMSP=GMSP)
              AUX[,SiteSN:=SSN[j]]
              AUX[,TSID:="ATH"]

              IM_SSN <- rbindlist(list(IM_SSN,AUX))

            }
          }
          OK <- !is.null(IM_SSN)
          if(OK){
            cat(sprintf("Done.\n"))

            IM_SSN[,RecordSN:=SET$RecordSN[i]]
            IM_SSN[,isRawRecord := FALSE]
            IM_SSN[,isRegularRecord := TRUE]
            IM_SSN[,TypeID := TID]
            IM_SSN[,TopLayerID := "outcrop"]
            # Write CSV ----
            cat(sprintf("> Writing S2O Intensity Table ..."))
            FILE <- file.path(IndexFolder,paste0(TID,".IntensityTable.",RUNID,".csv"))
            fwrite(IM_SSN,file=FILE,yaml = TRUE,append = TRUE)
            cat(sprintf(" Done!.\n"))
          }



        },
        error = function(e){
          cat(sprintf("* ERROR: %s *",e))
          OK <<- FALSE},
        warning = function(e){
          cat(sprintf("* WARNING: %s *",e))
          OK <<- FALSE}
      )

      if(!OK){cat(sprintf("* Skip record %s\n",SET$RecordSN[i]))}

    }
    else {
      cat(sprintf("* Skip (Record %s already processed for all sites)\n",SET$RecordSN[i]))}
  }

}



##
## -- 6. Database binding and Indexing (GMRI)---------------------------------------------------
# Copyright SRK (C) 2021
# Author: A.Verri K. (averri@srk.com.ar)
# Maintainer: A.Verri K. (averri@srk.com.ar)
# Date: 26-05-2021


.bindHeaderTables <- function(GMSP=NULL,OnlyKnownNEHRP=TRUE){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  ## Load Headers -------------------------------------------------

  PATH <- file.path(IndexFolder)


  cat(sprintf("> Binding HeaderTable \n "))
  PATTERN <- "*.HeaderTable.Rds"
  SET <- data.table(HeaderFilename=list.files(
    path = PATH,
    recursive = FALSE,
    pattern = PATTERN,
    full.names = FALSE,
    ignore.case=TRUE,
    include.dirs = FALSE))
  NR <- nrow(SET)
  HT <- data.table()
  stopifnot(NR>1)
  for(k in seq_len(NR)){
    cat(sprintf("> Binding header #%d/%d (Total Rows: %d)\n",k,NR,nrow(HT)))
    RDSfile <- file.path(PATH,SET$HeaderFilename[k])
    AUX <- readRDS(RDSfile)  |> as.data.table()
    HT <- rbindlist(list(HT,AUX),use.names=TRUE)
  }


  # Build TypeID --------------------------------------------------------
  HT[,TypeID:=str_sub(OwnerRecordSN,1,3)]

  # Tag Owner  -----------------------------------------------------------
  HT[ TypeID=="AT2", OwnerID:="NGA"]
  HT[ TypeID=="V2A", OwnerID:="GNS"]
  HT[ TypeID=="V2C", OwnerID:="CSN"]
  HT[ TypeID=="TRA", OwnerID:="GSC"]
  HT[ TypeID=="TRB", OwnerID:="CSN"]
  HT[ TypeID=="TRC", OwnerID:="ESM"]
  HT[ TypeID=="TRZ", OwnerID:="SRK"]
  HT[ TypeID=="LIS", OwnerID:="ICE"]
  HT[ TypeID=="ACB", OwnerID:="RED"]
  HT[ TypeID=="ACA", OwnerID:="IGP"]

  # Event data  -----------------------------------------------------------
  # PATH <- file.path(IndexFolder)
  # CSVfile <- file.path(PATH,"EventTable.csv")
  # stopifnot(file.exists(CSVfile))
  # ET <- read_csv(file=CSVfile,na=c("-"," ",""),  col_names = FALSE) |> as.data.table()
  # colnames(ET) <- ET[1] |> as.character()
  # ET <- ET[-1]
  # ET[,c("EventSN","")]

  # Station data  -----------------------------------------------------------

  FILE <- file.path(IndexFolder,"StationTable.Rds")
  stopifnot(file.exists(FILE))

  # ST <- fread(FILE,na.strings = c("--","-","")) #Check that StationTable.csv should not have NA VALUES
  ST <- readRDS(FILE)
  ST[Ts=="-" | Ts=="-1",Ts:=NA]
  ST[!is.na(Ts) & grepl(Ts,pattern = "<"),Ts:=str_remove(Ts,pattern = "<")]
  ST[!is.na(Ts) & grepl(Ts,pattern = ">"),Ts:=str_remove(Ts,pattern = ">")]
  ST[!is.na(Ts),Ts:=as.numeric(Ts)]
  ST[Vs30 =="-" | Vs30=="",Vs30 :=NA]
  ST[!is.na(Vs30),Vs30:=as.numeric(Vs30)]
  ST[NEHRP =="-" | NEHRP=="",NEHRP :=NA]
  ST[z1000=="-" | z1000=="",z1000:=NA]
  ST[!is.na(z1000),z1000:=as.numeric(z1000)]
  #  Fix NEHRP Ranges
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30<150,NEHRP:="E"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>150 & Vs30 <= 210,NEHRP:="DE"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>210 & Vs30 <= 300,NEHRP:="D"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>300 & Vs30 <= 440,NEHRP:="CD"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>440 & Vs30 <= 640,NEHRP:="C"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>640 & Vs30 <= 900,NEHRP:="BC"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>900 & Vs30 <= 1500,NEHRP:="B"]
  ST[!is.na(Vs30) & is.na(NEHRP) & Vs30>1500 ,NEHRP:="A"]


  # Tag LIS data  -----------------------------------------------------------
  HT.LIS <- HT[TypeID=="LIS"]
  cat(sprintf("> Removing %d unknown LIS Station Records\n",nrow(HT.LIS[is.na(OwnerStationSN)])))
  HT.LIS <- HT.LIS[!is.na(OwnerStationSN)]

  HT.LIS <- ST[TypeID=="LIS"][HT.LIS,on=c("OwnerStationSN","TypeID","OwnerID")]
  HT.LIS[is.na(StationLatitude) & !is.na(i.StationLatitude),StationLatitude:=i.StationLatitude]
  HT.LIS[is.na(StationLongitude) & !is.na(i.StationLongitude),StationLongitude:=i.StationLongitude]
  HT.LIS[is.na(StationName) & !is.na(i.StationName),StationName:=i.StationName]
  HT.LIS[,i.StationName:=NULL]
  HT.LIS[,i.StationLatitude:=NULL]
  HT.LIS[,i.StationLongitude:=NULL]

  cat(sprintf("> Parse Date & Time LIS \n"))
  # HT.LIS[is.na(EventTime),EventTime:="00:00:00"]
  # HT.LIS[,EventDate:=anytime(paste(EventDate, EventTime),asUTC = TRUE)]
  HT.LIS[,EventDate:=as_date(EventDate)]  |>  suppressWarnings()

  # Tag ACA data  -----------------------------------------------------------
  HT.ACA <- HT[TypeID=="ACA"]
  cat(sprintf("> Removing %d unknown ACA Station Records\n",nrow(HT.ACA[is.na(OwnerStationSN)])))
  HT.ACA <- HT.ACA[!is.na(OwnerStationSN)]

  HT.ACA <- ST[TypeID=="ACA"][HT.ACA,on=c("OwnerStationSN","TypeID","OwnerID")]
  HT.ACA[is.na(StationLatitude) & !is.na(i.StationLatitude),StationLatitude:=i.StationLatitude]
  HT.ACA[is.na(StationLongitude) & !is.na(i.StationLongitude),StationLongitude:=i.StationLongitude]
  HT.ACA[is.na(StationName) & !is.na(i.StationName),StationName:=i.StationName]
  HT.ACA[,i.StationName:=NULL]
  HT.ACA[,i.StationLatitude:=NULL]
  HT.ACA[,i.StationLongitude:=NULL]
  cat(sprintf("> Parse Date & Time ACA \n"))
  # HT.ACA[,EventDate:=anytime(EventDate,asUTC = TRUE)] # FIX Bug i PArseownerHeader
  HT.ACA[,EventDate:= as_date(EventDate)]  |>  suppressWarnings()
  # Tag ACB data  -----------------------------------------------------------
  HT.ACB <- HT[TypeID=="ACB"]
  cat(sprintf("> Removing %d unknown ACB Station Records\n",nrow(HT.ACB[is.na(OwnerStationSN)])))
  HT.ACB <- HT.ACB[!is.na(OwnerStationSN)]

  HT.ACB <- ST[TypeID=="ACB"][HT.ACB,on=c("OwnerStationSN","TypeID","OwnerID")]
  HT.ACB[is.na(StationLatitude) & !is.na(i.StationLatitude),StationLatitude:=i.StationLatitude]
  HT.ACB[is.na(StationLongitude) & !is.na(i.StationLongitude),StationLongitude:=i.StationLongitude]
  HT.ACB[is.na(StationName) & !is.na(i.StationName),StationName:=i.StationName]
  HT.ACB[,i.StationName:=NULL]
  HT.ACB[,i.StationLatitude:=NULL]
  HT.ACB[,i.StationLongitude:=NULL]

  cat(sprintf("> Parse Date & Time V2A \n"))
  # HT.ACB[grep(EventTime,pattern = "--"),EventTime:="00:00:00"]
  # HT.ACB[is.na(EventTime),EventTime:="00:00:00"]
  # HT.ACB[,EventDate:=anytime(paste(EventDate,EventTime),asUTC = TRUE)]
  HT.ACB[,EventDate:=as_date(EventDate)] |>  suppressWarnings()

  # Tag TRB data  -----------------------------------------------------------
  HT.TRB <- HT[TypeID=="TRB"]
  cat(sprintf("> Removing %d unknown TRA Station Records\n",nrow(HT.TRB[is.na(OwnerStationSN)])))
  HT.TRB <- HT.TRB[!is.na(OwnerStationSN)] # EventData only. No Station Data
  HT.TRB[,EventDate:=as_date(EventDate)] |>  suppressWarnings()

  # Tag TRA data  -----------------------------------------------------------
  HT.TRA <- HT[TypeID=="TRA"]
  cat(sprintf("> Removing %d unknown TRA Station Records\n",nrow(HT.TRA[is.na(OwnerStationSN)])))
  HT.TRA <- HT.TRA[!is.na(OwnerStationSN)] # EventData only. No Station Data
  HT.TRA[,EventDate:=as_date(EventDate)] |>  suppressWarnings()

  # Tag TRZ data  -----------------------------------------------------------
  HT.TRZ <- HT[TypeID=="TRZ"]
  cat(sprintf("> Removing %d unknown TRZ Station Records\n",nrow(HT.TRZ[is.na(OwnerStationSN)])))
  HT.TRZ <- HT.TRZ[!is.na(OwnerStationSN)] # EventData only. No Station Data
  HT.TRZ[,EventDate:= as_date(EventDate)] |>  suppressWarnings()
  # Tag TRC data  -----------------------------------------------------------
  HT.TRC <- HT[TypeID=="TRC"]
  cat(sprintf("> Removing %d unknown TRC Station Records\n",nrow(HT.TRC[is.na(OwnerStationSN)])))
  HT.TRC <- HT.TRC[!is.na(OwnerStationSN)]

  HT.TRC <- ST[TypeID=="TRC"][HT.TRC,on=c("OwnerStationSN","TypeID","OwnerID")]
  HT.TRC[is.na(StationLatitude) & !is.na(i.StationLatitude),StationLatitude:=i.StationLatitude]
  HT.TRC[is.na(StationLongitude) & !is.na(i.StationLongitude),StationLongitude:=i.StationLongitude]
  HT.TRC[is.na(StationName) & !is.na(i.StationName),StationName:=i.StationName]
  HT.TRC[,i.StationName:=NULL]
  HT.TRC[,i.StationLatitude:=NULL]
  HT.TRC[,i.StationLongitude:=NULL]

  cat(sprintf("> Parse Date & Time TRC \n"))
  # HT.TRC[is.na(EventTime),EventTime:="00:00:00"]
  # HT.TRC[,EventDate:=anytime(paste(EventDate, EventTime),asUTC = TRUE)]
  HT.TRC[,EventDate:= as_date(EventDate)] |>  suppressWarnings()


  # Tag V2A data  -----------------------------------------------------------
  HT.V2A <- HT[TypeID=="V2A"]
  cat(sprintf("> Removing %d unknown V2A Station Records\n",nrow(HT.V2A[is.na(OwnerStationSN)])))
  HT.V2A <- HT.V2A[!is.na(OwnerStationSN)]

  HT.V2A <- ST[TypeID=="V2A"][HT.V2A,on=c("OwnerStationSN","TypeID","OwnerID")]
  HT.V2A[is.na(StationLatitude) & !is.na(i.StationLatitude),StationLatitude:=i.StationLatitude]
  HT.V2A[is.na(StationLongitude) & !is.na(i.StationLongitude),StationLongitude:=i.StationLongitude]
  HT.V2A[is.na(StationName) & !is.na(i.StationName),StationName:=i.StationName]
  HT.V2A[,i.StationName:=NULL]
  HT.V2A[,i.StationLatitude:=NULL]
  HT.V2A[,i.StationLongitude:=NULL]

  cat(sprintf("> Parse Date & Time V2A \n"))
  # HT.V2A[is.na(EventTime),EventTime:="00:00:00"]
  # HT.V2A[,EventDate:=anytime(paste(EventDate,EventTime),asUTC = TRUE)]
  HT.V2A[,EventDate:= as_date(EventDate)] |>  suppressWarnings()


  # Tag AT2 Station data  -----------------------------------------------------------
  HT.AT2 <- HT[TypeID=="AT2"]

  HT.AT2[,OwnerStationSN:=NULL]
  FT <- .buildFileTable(TypeID="AT2",MULTI=FALSE)[,c("OwnerFileSN","RecordSN")]
  FT[,RSN:=sapply(OwnerFileSN,function(x){str_split(x,pattern = "[_]",simplify = TRUE)[1]})]
  FT[,OwnerFileSN:=NULL]
  idx <- duplicated(FT$RecordSN)
  FT <- FT[!idx]
  # FT <- unique(FT)

  PATH <- file.path(IndexFolder)
  FILE <- file.path(PATH,"AT2.RSN_SSN.csv")
  # OSN <- read_csv(FILE, col_names = FALSE) |> as.data.table()
  # colnames(OSN) <- OSN[1] |> as.character()
  # OSN <- OSN[-1]
  OSN <- fread(FILE,na.strings=c("-999"))[!is.na(OwnerStationSN)]
  OSN[,RSN:=paste0("RSN",RSN)]
  AUX <- OSN[FT,on="RSN"]
  AUX[,RSN:=NULL]
  HT.AT2 <- AUX[HT.AT2,on="RecordSN"]
  HT.AT2[,OwnerStationSN:=as.character(OwnerStationSN)]
  cat(sprintf("> Removing %d unknown AT2 Station Records\n",nrow(HT.AT2[is.na(OwnerStationSN)])))
  HT.AT2 <- HT.AT2[!is.na(OwnerStationSN)]

  HT.AT2 <- ST[TypeID=="AT2"][HT.AT2,on=c("OwnerStationSN","TypeID","OwnerID")]
  HT.AT2[is.na(StationLatitude) & !is.na(i.StationLatitude),StationLatitude:=i.StationLatitude]
  HT.AT2[is.na(StationLongitude) & !is.na(i.StationLongitude),StationLongitude:=i.StationLongitude]
  HT.AT2[is.na(StationName) & !is.na(i.StationName),StationName:=i.StationName]
  HT.AT2[,i.StationName:=NULL]
  HT.AT2[,i.StationLatitude:=NULL]
  HT.AT2[,i.StationLongitude:=NULL]
  # Tag AT2 Event data  -----------------------------------------------------------
  PATH <- file.path(IndexFolder)
  FILE <- file.path(PATH,"AT2.RSN_ESN.csv")
  stopifnot(file.exists(FILE))
  # OSN <- read_csv(FILE, col_names = FALSE) |> as.data.table()
  # colnames(OSN) <- OSN[1] |> as.character()
  # OSN <- OSN[-1]
  OSN <- fread(FILE,na.strings=c("-999",""))
  OSN <- OSN[!is.na(ESN)]
  OSN[,RSN:=paste0("RSN",RSN)]
  AUX <- OSN[FT,on="RSN"]
  AUX[,ESN:=NULL]

  cat(sprintf("> Parse Date & Time AT2 \n"))
  AUX[,YEAR:=str_pad(YEAR,side="left",width=4,pad="0")]
  AUX[,MODY:=str_pad(MODY,side="left",width=4,pad="0")]
  AUX[,HRMN:=str_pad(HRMN,side="left",width=4,pad="0")]
  AUX[!is.na(YEAR)& !is.na(MODY)& !is.na(HRMN),EventDate:=anytime(paste(YEAR,MODY,HRMN),asUTC = TRUE)]
  AUX[!is.na(YEAR)& !is.na(MODY)& is.na(HRMN),EventDate:=anytime(paste(YEAR,MODY),asUTC = TRUE)]
  AUX[!is.na(YEAR)& is.na(MODY) & is.na(HRMN),EventDate:=anytime(paste(YEAR),asUTC = TRUE)]
  AUX[,c("YEAR","MODY","HRMN"):=NULL]


  HT.AT2 <- AUX[HT.AT2,on="RecordSN"]
  HT.AT2[is.na(EpicenterLatitude) & !is.na(i.EpicenterLatitude),EpicenterLatitude:=i.EpicenterLatitude]
  HT.AT2[is.na(EpicenterLongitude) & !is.na(i.EpicenterLongitude),EpicenterLongitude:=i.EpicenterLongitude]
  HT.AT2[is.na(EpicenterDistance) & !is.na(i.EpicenterDistance),EpicenterDistance:=i.EpicenterDistance]
  HT.AT2[is.na(EventName) & !is.na(i.EventName),EventName:=i.EventName]
  HT.AT2[is.na(EventDate) & !is.na(i.EventDate) ,EventDate:=anytime(i.EventDate)|> as.Date()]
  HT.AT2[is.na(HypocenterDepth) & !is.na(i.HypocenterDepth),HypocenterDepth:=i.HypocenterDepth]
  HT.AT2[is.na(Magnitude) & !is.na(i.Magnitude),Magnitude:=i.Magnitude]
  HT.AT2[is.na(MagnitudeType) & !is.na(i.MagnitudeType),MagnitudeType:=i.MagnitudeType]
  HT.AT2[,(c("RSN","i.EpicenterLatitude","i.EpicenterLongitude","i.EpicenterDistance",
             "i.EventName","i.EventDate","i.HypocenterDepth","i.Magnitude","i.MagnitudeType")):=NULL]
  HT.AT2[,EventDate:= as_date(EventDate)] |>  suppressWarnings()


  rm(AUX,OSN,FT)
  # Bind Headers   --s---------------------------------------------------------
  HT <- rbindlist(list(HT.LIS,HT.ACA,HT.ACB,HT.AT2,HT.TRA,HT.TRB,HT.TRC,HT.TRZ,HT.V2A),use.names = TRUE,fill = TRUE)
  rm(HT.LIS,HT.ACA,HT.ACB,HT.AT2,HT.TRA,HT.TRB,HT.TRC,HT.TRZ,HT.V2A)

  # Clean EpicenterDistance   -----------------------------------------------------------
  DIGITS <- "[[:digit:]]+\\.*[[:digit:]]*"
  TEXT <- "[[:alpha:]]"
  HT[is.na(EpicenterDistance),EpicenterDistance:=-1]
  HT[grep(EpicenterDistance ,pattern = "km"),EpicenterDistance:=str_remove(EpicenterDistance,pattern = "km") ]
  HT[grep(EpicenterDistance,pattern =  TEXT),EpicenterDistance:=-1]
  GARB <- HT[is.na(as.numeric(EpicenterDistance)),EpicenterDistance] |> unique() |> suppressWarnings()
  HT[EpicenterDistance %in% GARB,EpicenterDistance:=-1]
  HT[,EpicenterDistance:=as.numeric(EpicenterDistance)]

  # Clean EpicenterLatitude  EpicenterLongitude ----
  HT[EpicenterLatitude=="--",EpicenterLatitude:=NA]
  HT[EpicenterLongitude=="--",EpicenterLongitude:=NA]

  # Clean HypocenterDepth   -----------------------------------------------------------
  HT[is.na(HypocenterDepth),HypocenterDepth:=-1]
  HT[ grep(HypocenterDepth ,pattern = "km"),HypocenterDepth :=str_remove(HypocenterDepth ,pattern = "km") ]
  HT[ grep(HypocenterDepth ,pattern = "R"),HypocenterDepth :=str_remove(HypocenterDepth ,pattern = "R") ]
  HT[grep(HypocenterDepth,pattern =  TEXT),HypocenterDepth:=-1]
  GARB <- HT[is.na(as.numeric(HypocenterDepth)),HypocenterDepth] |> unique() |> suppressWarnings()
  HT[HypocenterDepth %in% GARB,HypocenterDepth:=-1]
  HT[,HypocenterDepth:=as.numeric(HypocenterDepth)]

  # Clean Magnitude   -----------------------------------------------------------
  HT[is.na(Magnitude),Magnitude:=-1]
  HT[is.na(MagnitudeType),MagnitudeType:="U"]
  HT[ grep(Magnitude  ,pattern = "M"),Magnitude  :=str_remove(Magnitude  ,pattern = "M") ]
  GARB <- HT[is.na(as.numeric(Magnitude)),Magnitude]  |> unique()  |> suppressWarnings()
  HT[Magnitude %in% GARB,Magnitude:=-1]
  HT[,Magnitude:=as.numeric(Magnitude)|> round(digits = 1)]

  HT[grepl(MagnitudeType,pattern = DIGITS) & Magnitude==-1,`:=`(Magnitude=as.numeric(MagnitudeType) |> round(digits = 1))]
  HT[,MagnitudeType:=NULL]
  # Clean Ts   -----------------------------------------------------------
  HT[is.na(Ts),Ts:=-1]
  HT[Ts==">2",Ts:=runif(1,min=2,max=3)] # |> round(digits=1)]
  HT[Ts=="<0.1",Ts:=runif(1,min=0.05,max=0.1)]# |> round(digits=2)]
  GARB <- HT[is.na(as.numeric(Ts)),Ts] |> unique()
  HT[Ts %in% GARB,Ts:=-1]
  HT[,Ts:=as.numeric(Ts)]

  # Building Vs30 Ranges from Vs30 ------------------------
  HT[,EventTime:=NULL]
  ## Case 2: Vs30 == NA

  # Tag Records from NEHRP ---------------------

  cat(sprintf("> Building Vs30 Ranges from NEHRP\n"))
  HT[ NEHRP =="A",`:=`(VS30.min= 1500,VS30.max= 3000)]
  HT[ NEHRP =="B",`:=`(VS30.min= 900,VS30.max= 1500)]
  HT[ NEHRP =="BC",`:=`(VS30.min= 640,VS30.max= 900)]
  HT[ NEHRP =="C",`:=`(VS30.min= 440,VS30.max= 640)]
  HT[ NEHRP =="CD",`:=`(VS30.min= 300,VS30.max= 440)]
  HT[ NEHRP =="D",`:=`(VS30.min= 210,VS30.max= 300)]
  HT[ NEHRP =="DE",`:=`(VS30.min= 150,VS30.max= 210)]
  HT[ NEHRP =="E",`:=`(VS30.min= 10,VS30.max= 150)]

  cat(sprintf("> Building Vs30 Ranges from VS30 s\n"))
  HT[is.na(NEHRP)&  Vs30>= 1500 , `:=`(VS30.min= 1500,VS30.max= 3000,NEHRP="A")]
  HT[is.na(NEHRP)&  Vs30>= 900 & Vs30< 1500 , `:=`(VS30.min= 900,VS30.max= 1500,NEHRP="B")]
  HT[is.na(NEHRP)&  Vs30>= 640 & Vs30< 900 , `:=`(VS30.min= 640,VS30.max= 900,NEHRP="BC")]
  HT[is.na(NEHRP)&  Vs30>= 440 & Vs30< 640, `:=`(VS30.min= 440,VS30.max= 640,NEHRP="C")]
  HT[is.na(NEHRP)&  Vs30>= 300 & Vs30< 440, `:=`(VS30.min= 300,VS30.max= 440,NEHRP="CD")]
  HT[is.na(NEHRP)&  Vs30>= 210 & Vs30< 300, `:=`(VS30.min= 210,VS30.max= 300,NEHRP="D")]
  HT[is.na(NEHRP)&  Vs30>= 150 & Vs30< 210, `:=`(VS30.min= 150,VS30.max= 210,NEHRP="DE")]
  HT[is.na(NEHRP)&  Vs30< 150, `:=`(VS30.min= 10,VS30.max= 150,NEHRP="E")]




  # Tag Records from VS30 ---------------------
  cat(sprintf("> Tagging Vs30 ranges from CSCR data\n"))
  HT[ is.na(NEHRP) & CSCR =="S1",`:=`(VS30.min= 760,VS30.max= 1500,NEHRP="B")]
  HT[ is.na(NEHRP) & CSCR =="S2",`:=`(VS30.min= 350,VS30.max= 760,NEHRP="C")]
  HT[ is.na(NEHRP) & CSCR =="S3",`:=`(VS30.min= 180,VS30.max= 350,NEHRP="D")]
  HT[ is.na(NEHRP) & CSCR =="S4",`:=`(VS30.min= 60,VS30.max= 180,NEHRP="E")]


  cat(sprintf("> Tagging Vs30 ranges from EC8 data\n"))
  HT[ is.na(NEHRP) & EC8 =="A",`:=`(VS30.min= 800,VS30.max= 2500,NEHRP="B")]
  HT[ is.na(NEHRP) & EC8 =="B",`:=`(VS30.min= 360,VS30.max= 800,NEHRP="C")]
  HT[ is.na(NEHRP) & EC8 =="C",`:=`(VS30.min= 180,VS30.max=360,NEHRP="D")]
  HT[ is.na(NEHRP) & EC8 =="D",`:=`(VS30.min= 60,VS30.max= 180,NEHRP="E")]



  # TAg unknown regions
  HT[ TypeID=="V2A", `:=`(StationRegion="NEW ZEALAND",OwnerID="GNS")]
  HT[ TypeID=="TRA", `:=`(StationRegion="CANADA",OwnerID="GSC")]
  HT[ TypeID=="TRB", `:=`(StationRegion="CHILE",OwnerID="CSN")]
  HT[ TypeID=="TRZ", `:=`(StationRegion="AFRICA",OwnerID="SRK")]
  HT[ TypeID=="LIS", `:=`(StationRegion="CENTRAL AMERICA",OwnerID="ICE")]
  HT[ TypeID=="ACA", `:=`(StationRegion="PERU",OwnerID="IGP")]

  # *** REMOVE RECORDS WITH NO SITE DATA ****    --------------------------------------------------------

  # Check always that there are no records with other site data than NEHRP
  # ***************************************************************************
  if(OnlyKnownNEHRP){
    HT <- HT[!is.na(NEHRP)]
  }

  # ***************************************************************************

  # Rename Site Parameters   --------------------------------------------------------
  HT[,`:=`(GMX1=NULL,GMX3=NULL,SGS=NULL,NZS1170=NULL,NZS4203=NULL,EC8=NULL,CSCR=NULL,BSC=NULL,D2B=NULL)]
  # Remove EventName
  HT[,EventName:=NULL]
  setnames(HT,
           old=c("z1000","Vs30","NEHRP","Ts","VS30.min","VS30.max"),
           new=c("station.Z1000","station.VS30","station.NEHRP","station.Ts","station.VS30.min","station.VS30.max"))
  # Save --------------------------------------------------------
  cat(sprintf("> Save HeaderTable... "))
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"HeaderTable.Rds")
  saveRDS(HT,file=RDSfile,compress=TRUE)
  cat(sprintf("Done!\n"))
  return(HT[])
}

.bindOwnerIntensityTables <- function(){
  on.exit(expr={rm(list = ls())}, add = TRUE)

  IT <- data.table()
  # Load Raw Tables ---------------------------------------------------------------


  PATH <- file.path(IndexFolder)


  if(includeRawData==TRUE){
    PATTERN <- "*.OwnerIntensityTable.Rds"
    SET <- data.table(IMFilename=list.files(
      path = PATH,
      recursive = FALSE,
      pattern = PATTERN,
      full.names = TRUE,
      ignore.case=TRUE,
      include.dirs = FALSE))
    NR <- nrow(SET)
    cat(sprintf("> Binding raw data\n"))
    for(k in seq_len(NR)){
      cat(sprintf("> Binding table #%d/%d (Total Rows: %d)\n",k,NR,nrow(IT)))
      FILE <- SET$IMFilename[k]
      AUX <- readRDS(FILE)
      IT <- rbindlist(list(IT,AUX),use.names=TRUE,fill = TRUE)
    }
    IT[,`:=`(OwnerRecordSN=NULL)]
  }


  #
  # Load Regularized Tables ---------------------------------------------------------------
  PATTERN <- "*.RegularizedIntensityTable.*.csv"
  PATH <- file.path(IndexFolder)
  SET <- data.table(IMFilename=list.files(
    path = PATH,
    recursive = FALSE,
    pattern = PATTERN,
    full.names = TRUE,
    ignore.case=TRUE,
    include.dirs = FALSE))

  NR <- nrow(SET)
  cat(sprintf("> Binding regularized data \n"))
  for(k in seq_len(NR)){
    cat(sprintf("> Binding table #%d/%d (Total Rows: %d)\n",k,NR,nrow(IT)))
    FILE <- SET$IMFilename[k]
    AUX <- fread(FILE,yaml = TRUE)
    AUX[,`:=`(SiteSN="I")]
    IT <- rbindlist(list(IT,AUX),use.names=TRUE)
  }

  IT <- IT[TypeID!="TRA"]


  # Average Fs
  cat(sprintf("> Averaging regularized data by frequency\n"))
  IT <- IT[,.(
             TypeID,UN,
             IA=mean(IA),
             IA.cov=sd(IA)/mean(IA),

             IAu=mean(IAu),
             IAd=mean(IAd),

             PGA=mean(PGA),
             PGA.cov=sd(PGA)/mean(PGA),
             PGV=mean(PGV,na.rm=TRUE),
             PGD=mean(PGD,na.rm=TRUE),


             ARMS=mean(ARMS),
             ARMS.cov=sd( ARMS)/mean( ARMS),
             VRMS=mean(VRMS,na.rm=TRUE),

             DRMS=mean(DRMS,na.rm=TRUE),

             AZC=mean(AZC),
             VZC=mean(VZC,na.rm=TRUE),
             DZC=mean(DZC,na.rm=TRUE),

             ATo=mean(ATo),
             ATn=mean(ATn),
             VTo=mean(VTo),
             VTn=mean(VTn),
             DTo=mean(DTo),
             DTn=mean(DTn),

             CAV5=mean(CAV5),

             TmA=mean(TmA),

             DN05=mean(DN05),
             DN10=mean(DN10),
             DN15=mean(DN15),
             DN20=mean(DN20),
             DN25=mean(DN25),
             DN30=mean(DN30),
             DN40=mean(DN40),
             DN50=mean(DN50),
             DN75=mean(DN75),

             Dmax=mean(Dmax),
             D0595=mean(D0595),
             D2080=mean(D2080),
             D0575=mean(D0575),

             PPI=mean(PPI),
             EPI=mean(EPI),
             PDI=mean(PDI)


           ),by=.(RecordSN,OCID)] #all Fs


  # Check Integrity QA/QC ------------------------------------
  IT <- unique(IT,by=c("RecordSN","OCID" ))
  # Tag isUP  ------------------------------------------------
  cat(sprintf("> Tag Vertical Records (isUP) \n"))
  UP_ID <- c("UD","UP","V","HLZ","HNZ","DWN","DN","HHZ","Z","DOWN","VER","VERT","VERTICAL","VRT","BHZ","HHZ","HNZ","HLZ","HNZ","HLZ","HN3","HGZ","Z","UD","UP","V","Z","UPDO") |> unique()
  IT[toupper(OCID) %in% UP_ID, `:=`(isUP=TRUE,DIR="UP")]
  IT[!(toupper(OCID) %in% UP_ID), `:=`(isUP=FALSE,DIR="H2")]

  # Tag DIR  ------------------------------------------------
  cat(sprintf("> Tag H1 \n"))
  IT[DIR=="H2",isHmax:=(PGA==max(PGA)),by=c("RecordSN")] # all OCIDs
  IT[DIR=="H2" & isHmax==TRUE,DIR:="H1"]
  IT[,isHmax:=NULL]


  # TAG V/H ------------------------
  cat(sprintf("> Get V/H ratio... \n"))
  IT[order(DIR, decreasing = TRUE),`:=`(PGAV=PGA[1]),by=c("RecordSN")][,VH.PGA:=PGAV/PGA][,PGAV:=NULL] # all OCIDs
  IT[order(DIR, decreasing = TRUE),`:=`(ARMSV=ARMS[1]),by=c("RecordSN")][,VH.ARMS:=ARMSV/ARMS][,ARMSV:=NULL] # all OCIDs
  IT[order(DIR, decreasing = TRUE),`:=`(IAV=IA[1]),by=c("RecordSN")][,VH.IA:=IAV/IA][,IAV:=NULL] # all OCIDs






  # Load HeaderTable ----
  cat(sprintf("> Merge HeaderTable... "))
  # PATH <- file.path(IndexFolder,"Rds")
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"HeaderTable.Rds")
  HT <- readRDS(file=RDSfile)[,-c("StationName","station.Z1000","station.Ts","station.VS30.min","station.VS30.max")]
  setnames(HT,
           old=c("station.VS30","station.NEHRP"),
           new=c("OwnerVS30","OwnerNEHRP"))
  cat(sprintf("Done!\n"))
  COLS <- colnames(IT)[(colnames(IT) %in% colnames(HT))] #"RecordSN" "TypeID"
  IT <- HT[IT,on=COLS] |> unique(by=c("RecordSN","OCID"))
  rm(HT)

  IT[,`:=`(TopLayerID="station")]
  IT[,`:=`(SiteSN="I")]
  # Remove unknown sites ----------------------------------------
  IT <- IT[!is.na(OwnerVS30)][,OwnerVS30:=round(as.double(OwnerVS30))]
  # Save --------------------------------------------------------
  cat(sprintf("> Writting Itensity Table ..."))
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"OwnerIntensityTable.Rds")
  saveRDS(file=RDSfile,compress=TRUE, object=IT)
  # saveRDS(file=file.path(IndexFolder,"idx.Rds"),object=ITo[,.(RecordSN,SiteSN,Fs)] |> unique())
  cat(sprintf("Done!\n"))
  return(IT)
}

.bindIntensityTables <- function(GMSP=NULL){
  on.exit(expr={rm(list = ls())}, add = TRUE)


  # Load Tables ---------------------------------------------------------------

  # Get EXT
  PATTERN <- "*.OwnerIntensityTable.Rds"
  # Build Filetable
  cat(sprintf("> Indexing IntensityTables (raw data)\n"))
  # PATH <- file.path(IndexFolder,"Rds")
  PATH <- file.path(IndexFolder)

  SET <- data.table(IMFilename=list.files(
    path = PATH,
    recursive = FALSE,
    pattern = PATTERN,
    full.names = TRUE,
    ignore.case=TRUE,
    include.dirs = FALSE))



  NR <- nrow(SET)
  OIT <- data.table()
  if(NR>1){
    cat(sprintf("> Binding IntensityTables (raw data)\n"))
    for(k in seq_len(NR)){
      cat(sprintf("> Binding table #%d/%d (Total Rows: %d)\n",k,NR,nrow(OIT)))
      FILE <- SET$IMFilename[k]
      AUX <- readRDS(FILE)
      OIT <- rbindlist(list(OIT,AUX),use.names=TRUE)
    }
  }
  cat(sprintf("> Indexing IntensityTables (regularized data)\n"))
  # PATH <- file.path(IndexFolder,"Rds")
  PATH <- file.path(IndexFolder)
  PATTERN <- "*.IntensityTable.*.csv"
  SET <- data.table(IMFilename=list.files(
    path = PATH,
    recursive = FALSE,
    pattern = PATTERN,
    full.names = TRUE,
    ignore.case=TRUE,
    include.dirs = FALSE))

  list.files(
    path = PATH,
    recursive = FALSE,
    pattern = PATTERN,
    full.names = TRUE,
    ignore.case=TRUE,
    include.dirs = FALSE)

  NR <- nrow(SET)
  RIT <- data.table()
  if(NR>1){
    cat(sprintf("> Binding IntensityTables (regularized data)\n"))
    for(k in seq_len(NR)){
      cat(sprintf("> Binding table #%d/%d (Total Rows: %d)\n",k,NR,nrow(RIT)))
      FILE <- SET$IMFilename[k]
      AUX <- fread(FILE,yaml = TRUE)  |> suppressWarnings()
      # *** TEMPORARY FIX ***. FIX in owner intensity builder ***
      AUX[,EPI:=as.numeric(EPI)]
      RIT <- rbindlist(list(RIT,AUX),use.names=TRUE)
    }
  }
  cat(sprintf("> Binding IntensityTables\n"))
  # OIT[,`:=`(TSID=NULL,PGV=NULL,PGD=NULL,VZC=NULL,VRMS=NULL,VTo=NULL,VTn=NULL,DZC=NULL,DRMS=NULL,DTo=NULL,DTn=NULL)]
  # RIT[,`:=`(TSID=NULL,PGV=NULL,PGD=NULL,VZC=NULL,VRMS=NULL,VTo=NULL,VTn=NULL,DZC=NULL,DRMS=NULL,DTo=NULL,DTn=NULL)]
  #
  OIT[,`:=`(TSID=NULL,PGD=NULL,VRMS=NULL,DZC=NULL,DRMS=NULL,DTo=NULL,DTn=NULL)]
  RIT[,`:=`(TSID=NULL,PGD=NULL,VRMS=NULL,DZC=NULL,DRMS=NULL,DTo=NULL,DTn=NULL)]


  # *** TEMPORARY FIX ***. Remove TSID from builder
  # RIT[,TSID:=NULL]
  # OIT[,TSID:=NULL]
  # OIT[,BlockSN:=NULL]
  OIT[,OwnerRecordSN:=NULL]
  OIT[,TopLayerID:="station"]
  RIT[SiteSN=="I",TopLayerID:="station"]
  # Check Integrity QA/QC ------------------------------------
  # Remove owner records not processed in RIT
  RSN <- unique(RIT$RecordSN)
  OIT <- OIT[RecordSN %chin% RSN]
  OIT <- rbindlist(list(RIT,OIT),use.names=TRUE)
  rm(OIT,RIT)

  #IT <- IT[isRawRecord==FALSE | isRegularRecord==TRUE]
  IT <- unique(IT,by=c("RecordSN","OCID","SiteSN","TopLayerID","Fs"))

  # Tag isUP  ------------------------------------------------
  cat(sprintf("> Tag Vertical Records (isUP) \n"))
  UP_ID <- c("UD","UP","V","HLZ","HNZ","DWN","DN","HHZ","Z","DOWN","VER","VERT","VERTICAL","VRT","BHZ","HHZ","HNZ","HLZ","HNZ","HLZ","HN3","HGZ","Z","UD","UP","V","Z","UPDO") |> unique()
  IT[toupper(OCID) %in% UP_ID, `:=`(isUP=TRUE,DIR="UP")]
  IT[!(toupper(OCID) %in% UP_ID), `:=`(isUP=FALSE,DIR="H2")]

  # Tag DIR  ------------------------------------------------
  cat(sprintf("> Tag H1 \n"))
  IT[DIR=="H2",isHmax:=(PGA==max(PGA)),by=c("RecordSN","SiteSN","TopLayerID","isRawRecord","Fs")]
  IT[DIR=="H2" & isHmax==TRUE,DIR:="H1"]
  IT[,isHmax:=NULL]


  # TAG V/H ------------------------
  cat(sprintf("> Get V/H ratio... \n"))
  IT[order(DIR, decreasing = TRUE),`:=`(PGAV=PGA[1]),by=c("TopLayerID","RecordSN","SiteSN","isRawRecord","Fs")][,VH.PGA:=PGAV/PGA][,PGAV:=NULL]
  IT[order(DIR, decreasing = TRUE),`:=`(ARMSV=ARMS[1]),by=c("TopLayerID","RecordSN","SiteSN","isRawRecord","Fs")][,VH.ARMS:=ARMSV/ARMS][,ARMSV:=NULL]
  IT[order(DIR, decreasing = TRUE),`:=`(IAV=IA[1]),by=c("TopLayerID","RecordSN","SiteSN","isRawRecord","Fs")][,VH.IA:=IAV/IA][,IAV:=NULL]


  # Clean NA
  IT[is.na(CAV5),CAV5:=0]
  IT[is.na(EPI),CAV5:=0]
  IT[is.na(PPI),CAV5:=0]
  IT[is.na(PDI),CAV5:=0]

  IT[is.na(DN05),DN05:=0]
  IT[is.na(DN10),DN10:=0]
  IT[is.na(DN15),DN15:=0]
  IT[is.na(DN20),DN20:=0]
  IT[is.na(DN25),DN25:=0]
  IT[is.na(DN30),DN30:=0]
  IT[is.na(DN40),DN40:=0]
  IT[is.na(DN50),DN50:=0]
  IT[is.na(DN75),DN75:=0]

  # Get outcrop intensities ------------------------------------------------
  cat(sprintf("> Get outcrop intensities... \n"))
  AFo <- IT[TopLayerID=="outcrop",.(
    RecordSN,SiteSN,OCID,Fs=Fs,
    Dmax.outcrop=Dmax,
    D0595.outcrop=D0595,
    D2080.outcrop=D2080,
    VH.PGA.outcrop=VH.PGA,
    VH.IA.outcrop=VH.IA,
    VH.ARMS.outcrop=VH.ARMS,
    PGA.outcrop=PGA,
    PGV.outcrop=PGV,#PGD.outcrop=PGD,
    AZC.outcrop=AZC,
    VZC.outcrop=VZC,#DZC.outcrop=DZC,
    CAV5.outcrop=CAV5,
    IA.outcrop=IA,
    ARMS.outcrop=ARMS,
    # VRMS.outcrop=VRMS,DRMS.outcrop=DRMS,
    TmA.outcrop=TmA,
    PPI.outcrop=PPI,    EPI.outcrop=EPI,    PDI.outcrop=PDI,
    DN05.outcrop=DN05,DN10.outcrop=DN10,    DN15.outcrop=DN15,DN20.outcrop=DN20,
    DN25.outcrop=DN25,DN30.outcrop=DN30,    DN40.outcrop=DN40,DN50.outcrop=DN50,
    DN75.outcrop=DN75  )] |>  unique(by=c("RecordSN","OCID","SiteSN","Fs"))


  cat(sprintf("> Get station intensities... \n"))

  AFs <- IT[TopLayerID=="station" ,.(
    RecordSN,OCID,Fs,
    Dmax.station=Dmax,
    D0595.station=D0595,
    D2080.station=D2080,
    VH.PGA.station=VH.PGA,
    VH.IA.station=VH.IA,
    VH.ARMS.station=VH.ARMS,
    PGA.station=PGA,
    #PGD.station=PGD,
    PGV.station=PGV,
    AZC.station=AZC,
    VZC.station=VZC,#DZC.station=DZC,
    CAV5.station=CAV5,
    IA.station=IA,
    ARMS.station=ARMS,
    #VRMS.station=VRMS,#DRMS.station=DRMS,
    TmA.station=TmA,
    PPI.station=PPI,    EPI.station=EPI,    PDI.station=PDI,
    DN05.station=DN05,    DN10.station=DN10,    DN15.station=DN15,    DN20.station=DN20,
    DN25.station=DN25,DN30.station=DN30,DN40.station=DN40,DN50.station=DN50,
    DN75.station=DN75)]  |>   unique(by=c("RecordSN","OCID","Fs"))

  cat(sprintf("> Merge outcrop & bedrock intensities... \n"))
  AF <- AFs[AFo,on=c("RecordSN","OCID","Fs")]  |> na.omit()
  # AF[,i.SiteSN:=NULL]
  rm(AFo,AFs)


  # Load HeaderTable ----
  cat(sprintf("> Merge HeaderTable... "))
  # PATH <- file.path(IndexFolder,"Rds")
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"HeaderTable.Rds")
  # HT <- readRDS(file=RDSfile)[,-c("MagnitudeType","Magnitude","station.VS30.min",
  #                                 "EpicenterDistance","station.VS30.max","EventName",
  #                                 "StationName")]
  HT <- readRDS(file=RDSfile)[,-c("station.VS30.min","station.VS30.max","EventName","StationName")]
  setnames(HT,
           old=c("station.VS30","station.NEHRP","station.Z1000","station.Ts"),
           new=c("OwnerVS30","OwnerNEHRP","OwnerZ1000","OwnerTs"))
  cat(sprintf("Done!\n"))
  COLS <- colnames(AF)[(colnames(AF) %in% colnames(HT))]
  AF <- HT[AF,on=COLS]

  rm(HT)

  # Merge Site Data from Models ------------------------
  cat(sprintf("> Merge Site Data from Models... "))
  SiteTable <- readRDS(file=file.path(SitesFolder,"SiteTable.Rds"))
  setnames(SiteTable,old="SID",new="SiteSN")
  SiteTable[,D:=GMSP$D]
  COLS <- colnames(AF)[(colnames(AF) %in% colnames(SiteTable))]

  AF <- SiteTable[AF,on=COLS]  |> na.omit()
  colnames(IT) |> grep(pattern = "i\\.",value = TRUE)
  AF[,VS30:=as.double(VS30)]

  IT <-IT[,.(RecordSN,OCID,TypeID,DIR,UN,SFU,Scaled,NP,dt,Fs,isRawRecord,isRegularRecord,isUP)] |> unique(by=c("RecordSN","OCID","isRawRecord","Fs"))
  # rm(IT)
  COLS <- colnames(AF)[(colnames(AF) %in% colnames(IT))]
  ITo <- IT[AF,on=COLS]  |> na.omit()
  rm(AF,IT)
  colnames(ITo) |> grep(pattern = "i\\.",value = TRUE)
  cat(sprintf("Done!\n"))




  # Build amplification factors --------------------------------------

  cat(sprintf("> Build amplification factors... "))
  ITo[,`:=`(

    Ao2s.PGA=PGA.station/PGA.outcrop,
    Ao2s.IA=IA.station/IA.outcrop,
    Ao2s.CAV5=CAV5.station/CAV5.outcrop,
    Ao2s.ARMS=ARMS.station/ARMS.outcrop,
    Ao2s.AZC=AZC.station/AZC.outcrop,
    Ao2s.TmA=TmA.station/TmA.outcrop,
    Ao2s.PDI=PDI.station/PDI.outcrop,
    Ao2s.EPI=EPI.station/EPI.outcrop,
    Ao2s.PPI=PPI.station/PPI.outcrop

  )]


  # AF[,`:=`(
  #Ab2o.PGA=PGA.outcrop/PGA.bedrock,
  # Ab2s.PGA=PGA.station/PGA.bedrock,
  #   #Ab2o.PGV=PGV.outcrop/PGV.bedrock,
  #   # Ab2s.PGV=PGV.station/PGV.bedrock,
  #   # Ao2s.PGV=PGV.station/PGV.outcrop,
  #
  #   #Ab2o.PGD=PGD.outcrop/PGD.bedrock,
  #   # Ab2s.PGD=PGD.station/PGD.bedrock,
  #   # Ao2s.PGD=PGD.station/PGD.outcrop
  #   #Ab2o.IA=IA.outcrop/IA.bedrock,
  #   # Ab2s.IA=IA.station/IA.bedrock,
  #
  #
  #   #Ab2o.CAV5=CAV5.outcrop/CAV5.bedrock,
  #   # Ab2s.CAV5=CAV5.station/CAV5.bedrock,
  #
  # )]
  #
  # AF[,`:=`(
  #   #Ab2o.ARMS=ARMS.outcrop/ARMS.bedrock,
  #   # Ab2s.ARMS=ARMS.station/ARMS.bedrock,
  #
  #
  #   #Ab2o.VRMS=VRMS.outcrop/VRMS.bedrock,
  #   # Ab2s.VRMS=VRMS.station/VRMS.bedrock,
  #   # Ao2s.VRMS=VRMS.station/VRMS.outcrop,
  #
  #   #Ab2o.DRMS=DRMS.outcrop/DRMS.bedrock,
  #   # Ab2s.DRMS=DRMS.station/DRMS.bedrock,
  #   # Ao2s.DRMS=DRMS.station/DRMS.outcrop
  # )]
  #
  #
  # AF[,`:=`(
  #   #Ab2o.AZC=AZC.outcrop/AZC.bedrock,
  #   # Ab2s.AZC=AZC.station/AZC.bedrock,
  #
  #
  #   #Ab2o.TmA=TmA.outcrop/TmA.bedrock,
  #   # Ab2s.TmA=TmA.station/TmA.bedrock,
  #
  # )]
  #
  #
  # AF[,`:=`(
  #   #Ab2o.PDI=PDI.outcrop/PDI.bedrock,
  #   # Ab2s.PDI=PDI.station/PDI.bedrock,
  #
  #
  #   #Ab2o.EPI=EPI.outcrop/EPI.bedrock,
  #   # Ab2s.EPI=EPI.station/EPI.bedrock,
  #
  #
  #   #Ab2o.PPI=PPI.outcrop/PPI.bedrock,
  #   # Ab2s.PPI=PPI.station/PPI.bedrock,
  #
  #
  # )]

  setcolorder(ITo,neworder = c("RecordSN","SiteSN","NEHRP","Fs","VS30","VSa","Ts","Hs","HsR","HsF","HsC","Z1000","Z500","a","TypeID","OCID","DIR","isUP","isRegularRecord","isRawRecord","UN","SFU","Scaled","NP","dt"))

  cat(sprintf("Done!\n"))

  # Save --------------------------------------------------------
  cat(sprintf("> Writting Itensity Table ..."))
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"IntensityTable.Rds")
  system.time(saveRDS(file=RDSfile,compress=TRUE, object=ITo))
  saveRDS(file=file.path(IndexFolder,"idx.Rds"),object=ITo[,.(RecordSN,SiteSN,Fs)] |> unique())
  cat(sprintf("Done!\n"))

  return(ITo)
}

.buildIndexTable <- function(GMSP=NULL){
  # *************************************************************************
  # Existe un bug que produce NAs en OwnerID y otros campos que sí estan definidos en IT
  # Revisar codigo
  # *************************************************************************

  # Load IntensityTable ----
  cat(sprintf("> Load IntensityTable ... "))
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"IntensityTable.Rds")
  IT <- readRDS(file=RDSfile)
  COLS <- colnames(IT)[!(grepl(colnames(IT),pattern = "\\."))]
  IT <- IT[,..COLS]
  cat(sprintf("Done!\n"))
  #
  # AF <- AF[isRawRecord==FALSE,.(RecordSN,SiteSN,TypeID,OCID,DIR,isUP,isRegularRecord,StationName,StationRegion,StationLongitude,StationLatitude,EpicenterLongitude,EpicenterLatitude,EventDate  ,HypocenterDepth,Magnitude,MagnitudeType,EpicenterDistance)] |> unique()
  # AF[isRawRecord==TRUE,SiteSN:="I"]
  RecordSN_SET <- IT$RecordSN |> unique()
  # Load HeaderTable ----

  # cat(sprintf("> Load HeaderTable... "))
  # # PATH <- file.path(IndexFolder,"Rds")
  # PATH <- file.path(IndexFolder)
  # RDSfile <- file.path(PATH,"HeaderTable.Rds")
  # HT <- readRDS(file=RDSfile)
  # cat(sprintf("Done!\n"))

  # FileTable ----
  cat(sprintf("> Build FileTable... "))

  FT <- NULL
  ValidTypeID <- c("TRZ","TRA","TRB","TRC","V2A","ACA","ACB","LIS","AT2")
  for(TypeID in ValidTypeID){
    AUX <-.buildFileTable(TypeID)
    FT <- rbindlist(list(FT,AUX),use.names = TRUE)
  }
  rm(AUX)



  # remove dup columns and keep only records with inteNsities
  FT <- FT[RecordSN %chin%  unique(IT$RecordSN)]
  FT <- FT[,-c("OwnerStationSN")]
  # IT <- IT[,-c("OwnerRecordSN","OwnerStationSN","TypeID")]
  IT <- IT[,-c("OwnerStationSN","TypeID")]
  COLS <- colnames(FT)[colnames(FT) %in% colnames(IT)]

  MT1 <- IT[FT[is3D==FALSE],on=c("OwnerRecordSN", "RecordSN","OCID" )]
  MT2 <- IT[FT[is3D==TRUE],on=c("OwnerRecordSN","RecordSN" )]
  MT2[,i.OCID:=NULL]
  MT <- rbindlist(list(MT1,MT2),use.names = TRUE)
  rm(MT1,MT2,FT,IT)

   # Save --------------------------------------------------------
  cat(sprintf("> Save Data..."))
  PATH <- file.path(IndexFolder)
  RDSfile <- file.path(PATH,"MasterIndex.Rds")
  saveRDS(file=RDSfile,object=MT,compress=FALSE)
  cat(sprintf(".Done!\n"))
  return(MT)
}

.setNEHRPfromVS30 <- function(VS30){
  if(VS30>=1500){"A"}
  else if(VS30>=900 & VS30<1500){"B"}
  else if(VS30>=640 & VS30<900){"BC"}
  else if(VS30>=440 & VS30<640){"C"}
  else if(VS30>=300 & VS30<440){"CD"}
  else if(VS30>=210 & VS30<300){"D"}
  else if(VS30>=150 & VS30<210){"DE"}
  else if(VS30>=0 & VS30<150){"E"}
  else {NA}
}

.setEC8fromVS30 <- function(VS30){
  if(VS30<180){"D"}
  else  if(VS30>=180 & VS30<360){"C"}
  else  if(VS30>=360 & VS30<80){"B"}
  else  if(VS30>=800 ){"A"}
  else {NA}
}

.deg2dms <-  function (DEG, ROUND=FALSE){
  if(!is.na(DEG) & !is.nan(DEG) ){
    D <- as.integer(DEG)
    M_float <- 60*(DEG - D)
    M <- as.integer(M_float)
    S_float <- 60*(M_float - M)
    if(ROUND){
      S <- round(S_float, digits = 0 )
      ### After rounding, the seconds might become 60
      if (S == 60) {
        M <- M + 1
        S <- 0
      }
      if (M == 60){
        D <- D + 1
        M <- 0
      }
    } else {S <- S_float}
    DMS <- data.frame(D=D,M=M,S=S)
  } else {DMS <- data.frame(D=NA,M=NA,S=NA)}
  return (DMS)
}

.dms2deg <- function (DMS){
  DEG <- DMS$D + DMS$M/60+ DMS$S/3600
  return (DEG)
}

.getHEMI <- function(DMS,LAT=TRUE){
  # % Positive:'N' for longitudes and 'E' for latitudes
  # % Negative:'S' for longitudes and 'W' for latitudes
  if(DMS$D<0){
    HEMI <- ifelse(LAT,"S","W")
  } else {
    HEMI <- ifelse(LAT,"N","E")
  }
  return(HEMI)
}

.deg2SN <- function(DEG,LAT=TRUE){
  if(!is.na(DEG) & !is.nan(DEG) ){
    DMS <- deg2dms(DEG,ROUND=TRUE)
    NUM <- as.character(abs(DMS$S)+100*abs(DMS$M)+10000*abs(DMS$D)+1e7)
    n <- nchar(NUM)
    HEMI <- getHEMI(DMS,LAT)
    SN <- paste0(HEMI,substr(NUM,2,n))
    return(SN)
  }
  else return(NA)

}

.dms2SN <- function(DMS,LAT=TRUE){
  NUM <- as.character(abs(DMS$S)+100*abs(DMS$M)+10000*abs(DMS$D)+1e7)
  n <- nchar(NUM)
  HEMI <- getHEMI(DMS,LAT)
  SN <- paste0(HEMI,substr(NUM,2,n))
  return(SN)
}

.buildCoordinateSN <- function(LatitudeDEG,LongitudeDEG){
  if(!is.na(LatitudeDEG) & !is.na(LongitudeDEG)){
    return(paste0(deg2SN(LatitudeDEG,LAT=TRUE),   deg2SN(LongitudeDEG,LAT=FALSE)))
  }
  else{return(NA)}

}


## -- 7. Ground-Motion Export (GMEX)---------------------------------------------------
# Copyright SRK (C) 2021
# Author: A.Verri K. (averri@srk.com.ar)
# Maintainer: A.Verri K. (averri@srk.com.ar)
# Date: 18-06-2021
#

.taperRecord <- function(x){
  n <- length(x)
  cumIA <- cumsum((x*x)/as.numeric(x %*% x))
  iL_stop <- which(cumIA>=0.97)|> first()
  iL_pass <- which(cumIA>=0.95)|> first()
  if(iL_stop==iL_pass){
    iL_stop <- min(n,iL_pass+50)
  }
  LP <- .buildLowPassButtterworth(f=seq(1,n),Fstop = iL_stop, Fpass = iL_pass, Astop=0.01, Apass = 0.99)
  iH_stop <- which(cumIA>=0.03)|> first()
  iH_pass <- which(cumIA>=0.05)|> first()
  if(iH_stop==iH_pass){
    iH_stop <- max(1,iH_pass-50)
  }
  HP <- .buildHighPassButtterworth(f=seq(1,n),Fpass = iH_pass, Fstop = iH_stop, Astop=0.01, Apass = 0.99)
  return(HP*LP)
}

.exportRecords <- function(DT,RootPath=NULL,ExportPath=NULL){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  OK <- is.data.table(DT) && !is.null(RootPath) && !is.null(ExportPath)
  stopifnot(OK)

  ## xxxxxx -----------------------------------------------------------
  NR <- nrow(DT)
  for(k in seq_len(NR)){
    RSN <- DT$RecordSN[k]
    SSN <- DT$SiteSN[k]
    SF <- DT$SF[k]
    # RecordFilename <- FT[RecordSN== RSN]$RecordFilename
    RDSfile <- paste0(RootPath,DT$SourceFilename[k])
    GMSR <- readRDS(RDSfile)
    SiteSN <- names(GMSR$reg$TS)
    stopifnot(SSN %in% SiteSN)
    TS <- GMSR$reg$TS[[SSN]]
    TS[,W:=NULL]
    dt <- GMSR$reg$dt
    UN <- GMSR$reg$Units

    ## xxxxxx -----------------------------------------------------------
    COLS <- colnames(TS) |> grep(pattern = "^AT\\.",value = TRUE)
    AT <- TS[,..COLS]
    PGA <- apply(AT,2,function(x){max(abs(x))})
    W <- .taperRecord(AT[[which.max(PGA)]])
    TS[,colnames(TS):=lapply(.SD,function(x){x*W*SF})]

    ## Trim Zeros ---------------------------------------------------------------------
    TS <- unique(TS)
    NP <- nrow(TS)
    tn <- seq(from=0,to=NP-1)*dt
    TS <- cbind(tn,TS)
    TS <-rbind(TS,TS[1]*0)
    ## Save Record ---------------------------------------------------------------------

    CSVfile <- file.path(ExportPath,paste0(RSN,SSN,".csv"))
    RDSfile <- file.path(ExportPath,paste0(RSN,SSN,".Rds"))
    saveRDS(object=TS,compress=FALSE, file=RDSfile)
    write_csv(TS,CSVfile)
  }

}




## -- 5. Dynamic Site Response (GMSR) ---------------------------------------------------
# Build Site Profiles
# Copyright SRK (C) 2021
# Author & Maintainer: A.Verri K. (averri@srk.com.ar)
# Date: 26-06-2022

.buildGMSP <- function(
    Fs=100,NWo=2000,OVLP=75,D=0.05,Ao_stop=5,Ao_pass=50,TargetUnits="mm",
    Fpass_LPo = 20,  Fstop_LPo = 25,  Fpass_HPo = 0,  Fstop_HPo = 0){
  if(Fs>0){
    DownFs <-   Fs # Target frequency,
    UpFs <- min(300,4*Fs)#5*Fs
  }

  # PATH <- file.path(SitesFolder)
  # FilterTable <- readRDS(file=file.path(PATH,"FilterTable.Rds")) |> as.data.table()


  df <- Fs/NWo
  fs <- seq(from=0,by=df,length.out=NWo/2)

  Fpass_LP <- round(Fpass_LPo/df,digits=0)*df
  Fstop_LP <- round(Fstop_LPo/df,digits=0)*df
  Fpass_HP <- round(Fpass_HPo/df,digits=0)*df
  Fstop_HP <- round(Fstop_HPo/df,digits=0)*df

  # fs <- readRDS(file=FILE)#
  # stopifnot(NWo==2*length(fs))
  # if(SiteResponse==TRUE){
  #   PATH <- file.path(SitesFolder,"TFT")
  #   FILE <- file.path(PATH,paste0("fs.",Fs,".Rds"))
  #   stopifnot(file.exists(FILE))
  #   PATH <- file.path(SitesFolder,"TFT")
  #   FILE <- file.path(PATH,paste("S2O",Fs,NWo,"D",100*D,"Rds",sep="."))
  #   S2O <-  readRDS(file=FILE)
  #   FILE <- file.path(PATH,paste("S2B",Fs,NWo,"D",100*D,"Rds",sep="."))
  #   S2B <-  readRDS(file=FILE)
  # } else {
  #   S2O <- NULL
  #   S2B <- NULL
  # }


  # PATH <- file.path(SitesFolder)
  # SiteTable <- readRDS(file=file.path(PATH,"SiteTable.Rds")) |> as.data.table()
  #


  GMSP <-list(
    # SiteTable=SiteTable,
    TargetUnits = TargetUnits,
    # fs=fs,
    OVLP=OVLP ,
    NW=NWo,
    DownFs=DownFs,
    UpFs=UpFs,
    # df=df,
    # S2O=S2O,
    # S2B=S2B,
    # D=D,
    Ao_stop=Ao_stop,
    Ao_pass=Ao_pass,

    Fpass_HP = Fpass_HP,
    Fpass_LP = Fpass_LP,
    Fstop_HP = Fstop_HP,
    Fstop_LP = Fstop_LP

  )
  return(GMSP)
}

.setNEHRP <- function(Vs30){
  SC <- NULL
  if(Vs30>1500) {SC <- "A"}
  else if(Vs30>900 & Vs30<=1500){SC <- "B"}
  else if(Vs30>640 & Vs30<=900){SC <- "BC"}
  else if(Vs30>440 & Vs30<=640){SC <- "C"}
  else if(Vs30>300 & Vs30<=440){SC <- "CD"}
  else if(Vs30>210 & Vs30<=300){SC <- "D"}
  else if(Vs30>150 & Vs30<=210){SC <- "DE"}
  else SC <- "E"
  return(SC)
}

.getTs <-  function(SiteProfile) {
  DT <- SiteProfile
  NL <- nrow(DT)
  A <- vector(mode="double",length=NL-1)
  B <- vector(mode="double",length=NL-1)
  f <- vector(mode="double",length=NL)
  Hs <- sum(DT$hs)
  VSm <- DT$VSm
  hs <- DT$hs
  zm <- DT$zm
  for(j in seq(1,NL)){
    f[j+1] <- f[j]+hs[j]*(Hs-zm[j])/VSm[j]
    A[j] <- (VSm[j]*(f[j+1]-f[j]))^2/hs[j]
    B[j] <- hs[j]*(f[j+1]+f[j])^2
  }
  ws2 <- 4*sum(A)/sum(B)
  Ts <- 2*pi/sqrt(ws2)
  return(Ts)
}

.buildProfile <- function(Hs=NULL,Hr=NULL,mode="chaotic",SMP=NULL){
  zmax <- 0
  OV <- vector(mode="double",length=ceiling(2*(Hs+Hr)))
  OS <- vector(mode="character",length=ceiling(2*(Hs+Hr)))
  hs <- OV
  zo <- OV
  zi <- OV
  zm <- OV
  po <- OV
  pi <- OV
  pm <- OV
  FC <- OV
  UID <- OS
  WID <- OS
  GID <- OS

  MID <- OS
  FID <- OS
  OCR <- OV
  eo <- OV
  emax <- OV
  emin <- OV
  gs <- OV
  Gm <- OV
  VSm <- OV
  Dr <- OV

  # Build Soil Profiles
  # preconsolidation pressure
  pop <- runif(n=1,min=100,max=3000)

  NL <- 0
  while(sum(hs)<Hs){
    NL <- NL+1
    # hs[NL] <- min(runif(n=1,min=0.50,max=5.0),Hs-sum(hs)) |> round(digits=1)

    hs[NL] <- min(rtriangle(n=1,a=1.0,b=5.0,c=1.0),Hs-sum(hs)) |> round(digits=1)

    if(NL>1) {zo[NL] <- zi[NL-1]}
    zi[NL] <- zo[NL]+hs[NL]
    zm[NL] <- (zi[NL]+zo[NL])/2
  }
  # O <- vector(mode="double",length=NL+1)
  # S <- vector(mode="character",length=NL+1)
  Hs <- sum(hs)

  ## Soil Layers ----------------------------
  for(j in seq(1,NL)){
    SET <- c("F","S","G")
    if(mode=="chaotic" & j>10){
      # Set Material. Gravel,Fines or Sand?
      SET <- c(MID[1:NL][MID!=""],sample(SET,size=1))
    }

    # Set group ID
    MID[j] <- sample(SET,size=1)
    if(MID[j]=="G"){GID[j] <-"Gravels" }
    if(MID[j]=="S"){GID[j] <- "Sands"}
    if(MID[j]=="F"){GID[j] <- "Fines"}
    # Coarse Grained
    if(GID[j] %in% c("Gravels","Sands")){
      # Set Graduation
      SET <- c("W","P")
      if(mode=="chaotic" & j>10) {
        SET <- c(WID[1:NL][WID!=""],sample(SET,size=1))
      }
      WID[j] <- sample(SET,size=1)

      # Set Fines: Clayley mixes or silty mixes
      SET <-c("C","M")
      if(mode=="chaotic" & j>10) {
        SET <-  c(FID[1:NL][FID!="" & FID!="O"],sample(SET,size=1))
      }
      FID[j] <- sample(SET, size = 1)

      # Set FC
      FC[j] <- runif(n=1,min=5,max=50)

      if(FC[j]<5) {
        UID[j] <- paste0(MID[j],WID[j]) # GW,GP, SW,SP
        # Void Ratio
        NR <- VoidRatiosUSCS[USCS == UID[j],.N]
        stopifnot(NR==1)
        emin[j] <- VoidRatiosUSCS[USCS == UID[j],emin]
        emax[j] <- VoidRatiosUSCS[USCS == UID[j],emax]
      }
      if(FC[j]>=5 & FC[j]<=12)  {
        UID[j] <- paste0(MID[j],WID[j],"-",MID[j],FID[j]) # GW-GC
        # Void Ratio
        U1 <- paste0(MID[j],WID[j])
        NR <- VoidRatiosUSCS[USCS == U1,.N]
        stopifnot(NR==1)
        e1 <- VoidRatiosUSCS[USCS == U1,emin]
        e2 <- VoidRatiosUSCS[USCS == U1,emax]

        # Void Ratio
        U2 <- paste0(MID[j],FID[j])
        emin[j] <- min(e1,VoidRatiosUSCS[USCS == U2,emin])
        emax[j] <- max(e2,VoidRatiosUSCS[USCS == U2,emax])
      }
      if(FC[j]>12 ) {
        UID[j] <- paste0(MID[j],FID[j])
        # Void Ratio
        NR <- VoidRatiosUSCS[USCS == UID[j],.N]
        stopifnot(NR==1)
        emin[j] <- VoidRatiosUSCS[USCS == UID[j],emin]
        emax[j] <- VoidRatiosUSCS[USCS == UID[j],emax]
      }
      ## *************************************************************

      ## Void Ratios
      # eo[j] <- rtriangle(n=1,a=emin[j] ,b=emax[j] ,c=emin[j]+f/2*(emax[j]-emin[j]))
      #
      eo[j] <- runif(min=emin[j],max=emax[j],n=1)
      # Relative Density
      Dr[j] <- (emax[j]-eo[j])/(emax[j]-emin[j])

      # Unit Weight kN/m3
      gs[j] <- 17+5*Dr[j]

      # Overburden pressure (kPa)
      Ko <- 0.5
      if(j>1) {po[j] <- pi[j-1]}
      Dp <-  1/3*(1+2*Ko)*gs[j]*hs[j] #(* kPa*)
      pi[j] <- po[j] + Dp
      pm[j] <- (pi[j]+po[j])/2
      OCR[j] <- 1
      # Shear Module
      stopifnot(nrow(SMP[GroupID==GID[j]])>1)
      A <- SMP[GroupID==GID[j],A]
      Ce <- SMP[GroupID==GID[j],Ce]
      mo <- SMP[GroupID==GID[j],mo]
      Fe <- sapply(Ce,function(x){(x-eo[j])^2 / (1+eo[j])})
      G <- sapply(seq(1,length(Fe)),function(n){A[n]*Fe[n]*pm[j]^(mo[n])/1000})
      Gm[j] <- max(1,rnorm(n=1,mean=mean(G), sd=sd(G)) )
      VSm[j] <- sqrt(Gm[j]/gs[j])*100
    }

    # Fine Grained
    if(GID[j]=="Fines"){
      # Set Plasticity Flag and Clay/Silts
      SET <- c("O","C","M")
      if(mode=="chaotic" & j>9) {
        SET <-  c(FID[1:NL][FID!=""],sample(SET,size = 1))
      }
      FID[j] <- sample(SET, size = 1)

      if(zm[j]<=10){p <- 0.6}
      else if(zm[j]<=30){p <- 0.3}
      else {p <- 0.10}
      PID <- sample(c("H","L"),size=1, prob=c(p,1-p))
      # Set Liquid Limits
      if(PID=="H"){
        LL <- runif(n=1,min=50,max=100)
      } else {  #PID=="L"
        LL <- runif(n=1,min=8,max=50)
      }
      U_LINE <- 0.9*(LL-8) # >8 by def LL
      A_LINE <- ifelse(LL>=20,0.73*(LL-20),0)
      # Set Plasticity Index
      if(FID[j]=="C") {
        IP <- runif(n=1,min=A_LINE,max=U_LINE)
        gs[j] <- runif(n=1,min=16,max=20)}
      if(FID[j]=="M") {
        IP <- runif(n=1,min=0,max=A_LINE)
        gs[j] <- runif(n=1,min=16,max=20)}
      if(FID[j]=="O") {
        IP <- runif(n=1,min=0,max=A_LINE)
        gs[j] <- runif(n=1,min=15,max=16)
      }

      UID[j] <- paste0(FID[j],PID)

      # Set FC
      FC[j] <- runif(n=1,min=51,max=95)

      # Overburden pressure (kPa)
      Ko <- 0.5
      if(j>1) {po[j] <- pi[j-1]}
      Dp <-  1/3*(1+2*Ko)*gs[j]*hs[j] #(* kPa*)
      pi[j] <- po[j] + Dp
      pm[j] <- (pi[j]+po[j])/2

      OCR[j]=(pm[j]+pop)/pm[j]
      # Void Ratio
      NR <- VoidRatiosUSCS[USCS == UID[j],.N]
      stopifnot(NR==1)
      emin[j] <- VoidRatiosUSCS[USCS == UID[j],emin]
      emax[j] <- VoidRatiosUSCS[USCS == UID[j],emax]

      # Preconsolidation Pressure
      # eo[j] <- rtriangle(n=1,a=emin[j] ,b=emax[j] ,c=emin[j]+f/2*(emax[j]-emin[j]))
      if(OCR[j]>3){
        eo[j] <- rtriangle(n=1,a=emin[j],b=emax[j],c=emin[j])#runif(min=emin[j],max=emax[j],n=1)

      } else {
        eo[j] <- rtriangle(n=1,a=emin[j],b=emax[j])#runif(min=emin[j],max=emax[j],n=1)
      }
      # Relative Density
      # Dr[j] <- (emax[j]-eo[j])/(emax[j]-emin[j])
      # Shear Module
      stopifnot(nrow(SMP[GroupID==GID[j]])>1)

      A <- SMP[GroupID==GID[j],A]
      Ce <- SMP[GroupID==GID[j],Ce]
      mo <- SMP[GroupID==GID[j],mo]
      Fe <- sapply(Ce,function(x){(x-eo[j])^2 / (1+eo[j])})
      m <- approx(x=c(0,20,40,60,80,100),y=c(0,0.18,0.30,0.41,0.48,0.48),xout = IP)$y
      G <- sapply(seq(1,length(Fe)),function(n){A[n]*Fe[n]*(OCR[j]^m)*pm[j]^(mo[n])/1000})
      Gm[j] <- max(1,rnorm(n=1,mean=mean(G), sd=sd(G)) )
      VSm[j] <- sqrt(Gm[j]/gs[j])*100
    }
  }

  ## Rock layers -------------------------------------------------
  Hs <- sum(hs)
  while(sum(hs)<Hs+Hr){
    NL <- NL+1
    GID[NL] <- "Rock"
    hs[NL] <- 6# min(rtriangle(n=1,a=0.50,b=1,c=1.0),Hs+Hr-sum(hs)) |> round(digits=1)
    zo[NL] <- zi[NL-1]
    zi[NL] <- zo[NL]+hs[NL]
    zm[NL] <- (zi[NL]+zo[NL])/2
    gs[NL] <- runif(n=1,min=24,max=26)
    VSm[NL] <- runif(n=1,min=800,max=1500)
    Gm[NL] <- gs[NL] *(VSm[NL]/100)^2
  }

  ## Bedrock layer -------------------------------------------------
  gs[NL] <- 26
  VSm[NL] <- 2000
  Gm[NL] <- gs[NL] *(VSm[NL]/100)^2
  GID[NL] <- "Bedrock"


  PT <- data.table(
    USCS=UID[1:NL],
    GroupID=GID[1:NL],
    FC=FC[1:NL],
    zo=zo[1:NL],zi=zi[1:NL],zm=zm[1:NL],
    hs=hs[1:NL],
    eo=eo[1:NL],    emax=emax[1:NL],emin=emin[1:NL],
    OCR=OCR[1:NL],
    gs=gs[1:NL],
    pm=pm[1:NL],
    Gm= Gm[1:NL],
    VSm=VSm[1:NL])
  PT[GroupID!="Bedrock",dt:=hs/VSm]
  PT[GroupID=="Bedrock",dt:=0]
  PT[,VSa:=cumsum(hs)/cumsum(dt)]
  return(PT)
}

.buildProfileTable <- function(NSIM=1,Zmin=2,Zmax=120,Hr=30,NEHRP_Target=NULL,MODE="chaotic",SMP=ShearModelParameters){

  ProfileTable <- data.table()
  N <- 1
  I <- 1
  MAXIT <- 400
  CONTINUE <- N<NSIM & I<MAXIT
  while(CONTINUE){
    Hs <- runif(n=1,min=Zmin,max=Zmax) |> round(digits = 0)#sample(x=seq(from=ceiling(Zmin),to=ceiling(Zmax)),size = 1, replace = TRUE)

    PT <- .buildProfile(Hs=Hs,Hr=Hr,mode=MODE)
    VS30 <- PT[GroupID!="Bedrock" & zi<=30,.(VS30=30/sum(dt))]$VS30
    NEHRP <- .setNEHRP(VS30)
    I <- I+1
    OK <- is.null(NEHRP_Target) || (NEHRP==NEHRP_Target)
    if(OK){
      Hs <- PT[GroupID=="Bedrock",zi]
      VSa <- PT[,.(VSa=Hs/sum(dt))]$VSa
      Ts <- .getTs(PT[GroupID!="Bedrock"])
      AUX <- c(NEHRP,round(VSa,0),round(Hs,0),round(100*Ts))
      SID <- digest(object = AUX,algo="crc32")
      # Tag Profile
      PT[,SID:=SID]
      ProfileTable <- rbindlist(list(ProfileTable,PT))
      N <- N+1
      cat(sprintf("Found %d/%d NEHRP:%s Hs: %4.1f m VS30: %4.1f m/s VSa:%4.1f m/s \n",N,I,NEHRP,Hs,VS30,VSa))
    }
    CONTINUE <- (N<NSIM & I<MAXIT)
  }
  if(I>=MAXIT) cat(sprintf("Found %d sites until max iterations reached (MAXIT=%d)\n",N,I))
  return(ProfileTable)
}

.buildSiteTable <- function(ProfileTable){
  SiteTable <- data.table()
  IDX <- sample(unique(ProfileTable$SID))
  NS <- length(IDX)
  for(k in seq_len(NS)){
    SID <- IDX[k]
    # Get Profile properties
    PT <- ProfileTable[SID==IDX[k]][GroupID!="Bedrock"]
    # PT[,VSa:=cumsum(hs)/cumsum(dt)]
    Hs <- PT[,sum(hs)]
    HsR <- PT[GroupID=="Rock",sum(hs)] |> round(digits = 1)
    HsF <- PT[GroupID=="Fines",sum(hs)] |> round(digits = 1)
    HsC <- PT[GroupID %in% c("Gravels","Sands"),sum(hs)] #|> round(digits = 1)
    VSa <- PT[,last(VSa)]|> round(digits = 1)#PT[,sum(hs)/sum(dt)]#
    VS30 <- PT[zi<=30,30/sum(dt)]|> round(digits = 1)
    NEHRP <- .setNEHRP(VS30)
    Ts <- .getTs(PT)

    DATA <- PT[,.(zi,VSa)]
    MDL <- lm(data=DATA,formula = zi~VSa )
    Z500 <- predict.lm(MDL,newdata = data.frame(VSa=500))|> round(digits = 1)
    Z1000 <- predict.lm(MDL,newdata = data.frame(VSa=1000))|> round(digits = 1)

    a <-  Ts*VSa/Hs
    ST <- data.table(
      SID=SID,
      Hs=Hs,#round(Hs,digits = 1),
      HsR=HsR,
      HsF=HsF,
      HsC=HsC,
      Z1000=Z1000,
      Z500=Z500,
      VS30=VS30,#round(VS30,digits = 1),
      VSa=VSa,#round(VSm,digits = 1),
      Ts= Ts,#round(Ts,digits = 3),
      a = a,#round(a,digits = 3),
      NEHRP=NEHRP
    )
    SiteTable <- rbindlist(list(SiteTable,ST))
  }
  return(SiteTable)
}

.regularizeProfile <- function(ST=NULL,plot=FALSE){
  HIST <- hist(ST[Hs<=200,Hs],breaks = 200,plot = FALSE)
  stopifnot(all(HIST$counts)>0)
  NS <- floor(min(HIST$counts)/10)*10
  STR <- data.table()
  mids <- HIST$mids
  dHs <- 1/2*(mids[2]-mids[1])
  Hmin <- 0
  Hmax <- mids[1]-dHs
  STR <- rbindlist(list(STR,ST[Hs>Hmin & Hs<=Hmax][sample(NS)]))

  for(H in mids){
    Hmin <- H-dHs
    Hmax <- H+dHs
    STR <- rbindlist(list(STR,ST[Hs>Hmin & Hs<=Hmax][sample(NS)]))

  }
  if(plot==TRUE) hist(STR$Hs)
  return(STR)
}

.buildTransferFunction <- function(SiteTable,ProfileTable,FilterTable,Fs=125,NW=2000){
  df <- GMSP$df #Fs/NW #128 Hz: 0.0625  Hz
  f <-  seq(from=0,by=df,length.out=NW/2)
  NR <- nrow(SiteTable) # Number of Sites
  NF <- length(f) # Number of periods
  cat(sprintf("> Building %d Transfer Functions for %d periods...\n",NR,NF))
  ws  <- 2*pi*f
  # HT <- rbindlist(list(data.table(fs=f,H="A[NL,]",D=D),data.table(fs=f,H="B[NL,]",D=D)))
  HT <- data.table()

  for(k in seq_len(length.out = NR)){

    sid <- SiteTable$SID[k]
    cat(sprintf("Building Transfer function for site %s (%d/%d)\n",sid,k,NR))
    PT <- ProfileTable[SID==sid][GroupID!="Bedrock"]
    hs <- PT$hs
    rho <- PT$gs


    for(d in c(0.01,0.03,0.05,0.10,0.20)){
      NL <- nrow(PT)#(SiteTable$NL[k])-1
      A <- matrix(0, nrow=NL, ncol=NF)
      B <- matrix(0, nrow=NL, ncol = NF)
      A[1,] <- 1.0
      B[1,] <- 1.0

      VSi <- PT$VSm*sqrt(1+2*1i*d)
      for (i in seq_len(length.out = (NL-1))){
        alpha <- (rho[i]*VSi[i])/(rho[i+1]*VSi[i+1])
        Km  <- abs(ws)/VSi[i]
        A[i+1,] <- 0.5*(A[i,]*(1+alpha)*exp(1i*Km*hs[i])+B[i,]*(1-alpha)*exp(-1i*Km*hs[i]))
        B[i+1,] <- 0.5*(A[i,]*(1-alpha)*exp(1i*Km*hs[i])+B[i,]*(1+alpha)*exp(-1i*Km*hs[i]))
      }

      DT <-  data.table(A[NL,]) |> t() |> as.data.table()
      DT[,ID:=sid]
      DT[,m:="A[NL,]"]
      DT[,d:=d]
      HT <- rbindlist(list(HT,DT))
      DT <-  data.table(B[NL,]) |> t() |> as.data.table()
      DT[,ID:=sid]
      DT[,m:="B[NL,]"]
      DT[,d:=d]
      HT <- rbindlist(list(HT,DT))
      rm(DT,A,B)
    } # end for d
  } # end for k
  # DT <-  data.table(f) |> t() |> as.data.table()
  # DT[,ID:="00000000"]
  # DT[,m:="fs"]
  # DT[,d:=0]
  # HT <- rbindlist(list(HT,DT))
  setcolorder(HT,c("ID","m","d"))
  return(HT)
}

.buildHomogeneousProfile <- function(Hs=NULL, pop= runif(n=1,min=100,max=3000)){
  zmax <- 0
  OV <- vector(mode="double",length=ceiling(2*(Hs+Hr)))
  OS <- vector(mode="character",length=ceiling(2*(Hs+Hr)))
  hs <- OV
  zo <- OV
  zi <- OV
  zm <- OV
  po <- OV
  pi <- OV
  pm <- OV
  FC <- OV
  UID <- OS
  WID <- OS
  GID <- OS

  MID <- OS
  FID <- OS
  OCR <- OV
  eo <- OV
  emax <- OV
  emin <- OV
  gs <- OV
  Gm <- OV
  VSm <- OV
  Dr <- OV

  # Build Soil Profiles
  SMP <- ShearModelParameters
  NL <- 0

  while(sum(hs)<Hs){
    NL <- NL+1
    hs[NL] <- Hs/NL |> round(digits=1)
    if(NL>1) {zo[NL] <- zi[NL-1]}
    zi[NL] <- zo[NL]+hs[NL]
    zm[NL] <- (zi[NL]+zo[NL])/2
  }

  Hs <- sum(hs)

  ## Soil Layers ----------------------------
  for(j in seq(1,NL)){
    SET <- c("F","S","G")
    if(mode=="chaotic" & j>10){
      # Set Material. Gravel,Fines or Sand?
      SET <- c(MID[1:NL][MID!=""],sample(SET,size=1))
    }

    # Set group ID
    MID[j] <- sample(SET,size=1)
    if(MID[j]=="G"){GID[j] <-"Gravels" }
    if(MID[j]=="S"){GID[j] <- "Sands"}
    if(MID[j]=="F"){GID[j] <- "Fines"}
    # Coarse Grained
    if(GID[j] %in% c("Gravels","Sands")){
      # Set Graduation
      SET <- c("W","P")
      if(mode=="chaotic" & j>10) {
        SET <- c(WID[1:NL][WID!=""],sample(SET,size=1))
      }
      WID[j] <- sample(SET,size=1)

      # Set Fines: Clayley mixes or silty mixes
      SET <-c("C","M")
      if(mode=="chaotic" & j>10) {
        SET <-  c(FID[1:NL][FID!="" & FID!="O"],sample(SET,size=1))
      }
      FID[j] <- sample(SET, size = 1)

      # Set FC
      FC[j] <- runif(n=1,min=5,max=50)

      if(FC[j]<5) {
        UID[j] <- paste0(MID[j],WID[j]) # GW,GP, SW,SP
        # Void Ratio
        NR <- VoidRatiosUSCS[USCS == UID[j],.N]
        stopifnot(NR==1)
        emin[j] <- VoidRatiosUSCS[USCS == UID[j],emin]
        emax[j] <- VoidRatiosUSCS[USCS == UID[j],emax]
      }
      if(FC[j]>=5 & FC[j]<=12)  {
        UID[j] <- paste0(MID[j],WID[j],"-",MID[j],FID[j]) # GW-GC
        # Void Ratio
        U1 <- paste0(MID[j],WID[j])
        NR <- VoidRatiosUSCS[USCS == U1,.N]
        stopifnot(NR==1)
        e1 <- VoidRatiosUSCS[USCS == U1,emin]
        e2 <- VoidRatiosUSCS[USCS == U1,emax]

        # Void Ratio
        U2 <- paste0(MID[j],FID[j])
        emin[j] <- min(e1,VoidRatiosUSCS[USCS == U2,emin])
        emax[j] <- max(e2,VoidRatiosUSCS[USCS == U2,emax])
      }
      if(FC[j]>12 ) {
        UID[j] <- paste0(MID[j],FID[j])
        # Void Ratio
        NR <- VoidRatiosUSCS[USCS == UID[j],.N]
        stopifnot(NR==1)
        emin[j] <- VoidRatiosUSCS[USCS == UID[j],emin]
        emax[j] <- VoidRatiosUSCS[USCS == UID[j],emax]
      }
      ## *************************************************************

      ## Void Ratios
      # eo[j] <- rtriangle(n=1,a=emin[j] ,b=emax[j] ,c=emin[j]+f/2*(emax[j]-emin[j]))
      #
      eo[j] <- runif(min=emin[j],max=emax[j],n=1)
      # Relative Density
      Dr[j] <- (emax[j]-eo[j])/(emax[j]-emin[j])

      # Unit Weight kN/m3
      gs[j] <- 17+5*Dr[j]

      # Overburden pressure (kPa)
      Ko <- 0.5
      if(j>1) {po[j] <- pi[j-1]}
      Dp <-  1/3*(1+2*Ko)*gs[j]*hs[j] #(* kPa*)
      pi[j] <- po[j] + Dp
      pm[j] <- (pi[j]+po[j])/2
      OCR[j] <- 1
      # Shear Module
      stopifnot(nrow(SMP[GroupID==GID[j]])>1)
      A <- SMP[GroupID==GID[j],A]
      Ce <- SMP[GroupID==GID[j],Ce]
      mo <- SMP[GroupID==GID[j],mo]
      Fe <- sapply(Ce,function(x){(x-eo[j])^2 / (1+eo[j])})
      G <- sapply(seq(1,length(Fe)),function(n){A[n]*Fe[n]*pm[j]^(mo[n])/1000})
      Gm[j] <- max(1,rnorm(n=1,mean=mean(G), sd=sd(G)) )
      VSm[j] <- sqrt(Gm[j]/gs[j])*100
    }

    # Fine Grained
    if(GID[j]=="Fines"){
      # Set Plasticity Flag and Clay/Silts
      SET <- c("O","C","M")
      if(mode=="chaotic" & j>9) {
        SET <-  c(FID[1:NL][FID!=""],sample(SET,size = 1))
      }
      FID[j] <- sample(SET, size = 1)

      if(zm[j]<=10){p <- 0.6}
      else if(zm[j]<=30){p <- 0.3}
      else {p <- 0.10}
      PID <- sample(c("H","L"),size=1, prob=c(p,1-p))
      # Set Liquid Limits
      if(PID=="H"){
        LL <- runif(n=1,min=50,max=100)
      } else {  #PID=="L"
        LL <- runif(n=1,min=8,max=50)
      }
      U_LINE <- 0.9*(LL-8) # >8 by def LL
      A_LINE <- ifelse(LL>=20,0.73*(LL-20),0)
      # Set Plasticity Index
      if(FID[j]=="C") {
        IP <- runif(n=1,min=A_LINE,max=U_LINE)
        gs[j] <- runif(n=1,min=16,max=20)}
      if(FID[j]=="M") {
        IP <- runif(n=1,min=0,max=A_LINE)
        gs[j] <- runif(n=1,min=16,max=20)}
      if(FID[j]=="O") {
        IP <- runif(n=1,min=0,max=A_LINE)
        gs[j] <- runif(n=1,min=15,max=16)
      }

      UID[j] <- paste0(FID[j],PID)

      # Set FC
      FC[j] <- runif(n=1,min=51,max=95)

      # Overburden pressure (kPa)
      Ko <- 0.5
      if(j>1) {po[j] <- pi[j-1]}
      Dp <-  1/3*(1+2*Ko)*gs[j]*hs[j] #(* kPa*)
      pi[j] <- po[j] + Dp
      pm[j] <- (pi[j]+po[j])/2

      OCR[j]=(pm[j]+pop)/pm[j]
      # Void Ratio
      NR <- VoidRatiosUSCS[USCS == UID[j],.N]
      stopifnot(NR==1)
      emin[j] <- VoidRatiosUSCS[USCS == UID[j],emin]
      emax[j] <- VoidRatiosUSCS[USCS == UID[j],emax]

      # Preconsolidation Pressure
      # eo[j] <- rtriangle(n=1,a=emin[j] ,b=emax[j] ,c=emin[j]+f/2*(emax[j]-emin[j]))
      if(OCR[j]>3){
        eo[j] <- rtriangle(n=1,a=emin[j],b=emax[j],c=emin[j])#runif(min=emin[j],max=emax[j],n=1)

      } else {
        eo[j] <- rtriangle(n=1,a=emin[j],b=emax[j])#runif(min=emin[j],max=emax[j],n=1)
      }
      # Relative Density
      # Dr[j] <- (emax[j]-eo[j])/(emax[j]-emin[j])
      # Shear Module
      stopifnot(nrow(SMP[GroupID==GID[j]])>1)

      A <- SMP[GroupID==GID[j],A]
      Ce <- SMP[GroupID==GID[j],Ce]
      mo <- SMP[GroupID==GID[j],mo]
      Fe <- sapply(Ce,function(x){(x-eo[j])^2 / (1+eo[j])})
      m <- approx(x=c(0,20,40,60,80,100),y=c(0,0.18,0.30,0.41,0.48,0.48),xout = IP)$y
      G <- sapply(seq(1,length(Fe)),function(n){A[n]*Fe[n]*(OCR[j]^m)*pm[j]^(mo[n])/1000})
      Gm[j] <- max(1,rnorm(n=1,mean=mean(G), sd=sd(G)) )
      VSm[j] <- sqrt(Gm[j]/gs[j])*100
    }
  }

  ## Rock layers -------------------------------------------------
  Hs <- sum(hs)
  while(sum(hs)<Hs+Hr){
    NL <- NL+1
    GID[NL] <- "Rock"
    hs[NL] <- 6# min(rtriangle(n=1,a=0.50,b=1,c=1.0),Hs+Hr-sum(hs)) |> round(digits=1)
    zo[NL] <- zi[NL-1]
    zi[NL] <- zo[NL]+hs[NL]
    zm[NL] <- (zi[NL]+zo[NL])/2
    gs[NL] <- runif(n=1,min=24,max=26)
    VSm[NL] <- runif(n=1,min=800,max=1500)
    Gm[NL] <- gs[NL] *(VSm[NL]/100)^2
  }

  ## Bedrock layer -------------------------------------------------
  gs[NL] <- 26
  VSm[NL] <- 2000
  Gm[NL] <- gs[NL] *(VSm[NL]/100)^2
  GID[NL] <- "Bedrock"


  PT <- data.table(
    USCS=UID[1:NL],
    GroupID=GID[1:NL],
    FC=FC[1:NL],
    zo=zo[1:NL],zi=zi[1:NL],zm=zm[1:NL],
    hs=hs[1:NL],
    eo=eo[1:NL],    emax=emax[1:NL],emin=emin[1:NL],
    OCR=OCR[1:NL],
    gs=gs[1:NL],
    pm=pm[1:NL],
    Gm= Gm[1:NL],
    VSm=VSm[1:NL])
  PT[GroupID!="Bedrock",dt:=hs/VSm]
  PT[GroupID=="Bedrock",dt:=0]
  PT[,VSa:=cumsum(hs)/cumsum(dt)]
  return(PT)
}

.predictKh <- function(.SD){
  stopifnot(nrow(.SD)==1)

  OwnerNEHRP <- .SD$OwnerNEHRP
  SC <- .SD$SC
  PGA <- .SD$PGA
  PGV <- .SD$PGV
  ARMS <- .SD$ARMS
  IA <- .SD$IA
  CAV5 <- .SD$CAV5
  Dmax <- .SD$Dmax
  PPI <- .SD$PPI
  EPI <- .SD$EPI
  PDI <- .SD$PDI
  D0595 <- .SD$D0595
  AZC <- .SD$AZC
  TypeID <- .SD$TypeID
  UN <- .SD$UN
  isUP <- .SD$isUP
  VH.PGA <- .SD$VH.PGA
  Magnitude <- .SD$Magnitude
  Dn <- c(.SD$DN05,.SD$DN10,.SD$DN15,.SD$DN20,.SD$DN25,.SD$DN30,.SD$DN40,.SD$DN50,.SD$DN75)
  ky <- c(0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50,0.75)*PGA
  DATA <- data.table(Dn,ky)
  DATA <- DATA[Dn>0][,.(LnKy=log(ky),LnDn=log(Dn),Dn)]
  MDL <- lm(data=DATA,LnKy~.)
  R2a <- summary(MDL)$adj.r.squared |> round(digits = 3)
  Da <- c(5,10,15,20,25,30,45,60,75,100,125,150,175,200,250,300,350,400)
  kmax <- predict(MDL,newdata = data.table(Dn=Da,LnDn=log(Da))) |> unname() |> exp()
  kh <- kmax/PGA
  return(data.table(OwnerNEHRP,SC,PGA,PGV,ARMS,IA,CAV5,Dmax,PPI,EPI,PDI,D0595,AZC,TypeID,UN,isUP,Da,VH.PGA,Magnitude,kmax=round(kmax,digits=4),kh=round(kh,digits=3),R2a))
}
