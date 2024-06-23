
#' Title
#'
#' @param .x data.table
#' @param s vector
#' @param ts vector
#' @param Units character
#' @param kh numeric
#'
#' @return
#' @export
#'
#' @examples
get_nd <- function(.x=NULL,s=NULL,ts=NULL,Units="mm",kh){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  # check null .x and s
  if(is.null(s) & is.null(.x)){
    stop("No acceleration data provided.")
  }
  # check NULLs
  if(is.null(ts) & is.null(.x)){
    stop("No time data provided.")
  }
  g <- .getG(Units)
  if(!is.null(.x)){
    ts <- .x$t
    s <- .x$s
  }
  dt <- mean(diff(ts))
  nd <- .getDN(s=s,dt=dt,kh=kh, g=g) 
  return(nd)
}

