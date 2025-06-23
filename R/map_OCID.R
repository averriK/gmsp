#' Title
#'
#' @param .x data.table
#'
#' @return list
#' @export
#'
#' @import data.table

map_OCID <- function(.x) {
    on.exit(expr = {
        rm(list = ls())
    }, add = TRUE)
    . <- NULL
    DT <- copy(.x)
    stopifnot(names(DT) %in% c("OCID", "AT"))
    UP_ID <- c("BHZ", "HHZ", "HLZ", "HNZ", "HN3", "HGZ", "DWN", "DN", "Z", "DOWN", "V", "VER", "VERT", "VERTICAL", "VRT", "UD", "UP", "UPDO") |> unique()
    OD <- DT[, .(PGA = max(abs(AT)), isUP = (toupper(OCID) %in% toupper(UP_ID)), DIR = ""), by = .(OCID)]
    OD[isUP == TRUE, DIR := "UP"]
    OD[isUP == FALSE, DIR := ifelse(PGA == max(PGA), "H1", "H2")]
    OD[, `:=`(isUP = NULL, PGA = NULL)]

    return(OD)
}
