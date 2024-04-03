devtools::load_all()
library(data.table)
RecordsFolder <- file.path("/Users/averri/Database/gmdb/source/tables")
FILES <- list.files(RecordsFolder,pattern = "*.Rds",full.names = FALSE)
for(FILE in FILES){
  SET <- readRDS(file.path(RecordsFolder,FILE))
  SET <- lapply(SET, function(x) {
    x$Header <- NULL
    return(x)

  })
  BASE <- tools::file_path_sans_ext(basename(FILE))
  FILE <- file.path(RecordsFolder,paste0(BASE,"_R.Rds"))
  saveRDS(SET,file=FILE)
#
}
