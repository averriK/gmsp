#' Title
#'
#' @param .x data.table
#'
#' @return
#' @export
#' 
#' @import data.table 
#'
#' @examples
map_OCID <- function(.x){
  DT <- copy(.x)
  stopifnot(names(DT) %in% c("OCID","AT"))
  OD <- DT[,.(PGA=max(abs(AT))),by=.(OCID)]
  UP_ID <- c("BHZ","HHZ","HLZ","HNZ","HN3","HGZ","DWN","DN","Z","DOWN","V","VER","VERT","VERTICAL","VRT","UD","UP","UPDO") |> unique()
  OD[,isUP:=(OCID %in% UP_ID)[1]]
  
  OD[,DIR:="H2"]
  OD[isUP==TRUE,DIR:="UP"]
  OD[!isUP & PGA==max(PGA),DIR:="H1"]
  OD[!isUP & PGA==min(PGA),DIR:="H2"]
  OD[,`:=`(isUP=NULL,PGA=NULL)]
  
  return(OD)
}
