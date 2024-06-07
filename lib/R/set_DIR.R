#' Title
#'
#' @param AT data.table
#'
#' @return
#' @export
#' 
#' @importFrom data.table function
#'
#' @examples
set_DIR <- function(AT){
  OCID <- colnames(AT) |> unique()
  if(length(OCID)!=3){
    stop("Invalid record. No vertical records") # invalid record. No vertical records
  }
  DIR <- OCID
  UP_ID <- c("BHZ","HHZ","HLZ","HNZ","HN3","HGZ","DWN","DN","Z","DOWN","V","VER","VERT","VERTICAL","VRT","UD","UP","UPDO") |> unique()
  idx_UP <- toupper(OCID)%in%UP_ID
  if(sum(idx_UP)!=1){
    stop("Unknown UP tag") # invalid record. No vertical records or multiple matches for UP_ID
  }
  DIR[idx_UP] <- "UP"
  PGA <- AT[,sapply(.SD,function(x){max(abs(x))})]
  idx_H1 <- PGA==max(PGA[!idx_UP])
  if(sum(idx_H1)==1) {
    DIR[idx_H1] <- "H1"
    DIR[!idx_H1 & !idx_UP] <- "H2"
  } else {
    DIR[idx_H1[1]] <- "H1"
    DIR[idx_H1[2]] <- "H2"
  }
  # setnames(AT,old=OCID,new=DIR)
  return(data.table(OCID=OCID,DIR=DIR))
}
