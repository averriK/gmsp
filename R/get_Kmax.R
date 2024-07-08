#' Title
#'
#' @param .SD data.table
#' @param TargetUnits character
#'
#' @return
#' @export
#'
#' @examples
get_kmax <- function(.SD,TargetUnits="mm",interval=TRUE){
  AUX <- .SD[,.(PGA,ky=as.numeric(gsub("DN", "", ID))*PGA/100,Dn=round(value,4))][Dn>0 & ky>0]
  DATA <- AUX[,.(LnKy=log(ky),LnDn=log(Dn),Dn,IDn=1/Dn)]
  MDL <- lm(data=DATA,LnKy~.)
  
  Da <- c(1,2.5,5,7.5,10,12.5,15,20,25,30,35,40,50,60,70,80,90,100,125,150,175,200,250,300,350,400,500,750,1000,1500,2000,2500,3000,3500,4000,5000,7500,10000)
  LnKmax <- predict(MDL,newdata = data.table(LnDn=log(Da),Dn=Da,IDn=1/Da),interval="prediction") 
  
  
  kmax <- exp(LnKmax[,"fit"])  |> round(digits = 3)
  kh <- (100*kmax/PGA) |> round(digits = 1)
  # Model performance
  Yp <- predict(MDL) 
  RSS <- DATA$LnKy - Yp
  RMSE <- sqrt(mean(RSS^2))
  R2 <- summary(MDL)$adj.r.squared |> round(digits = 3)
  if(interval==TRUE){
    kmax95 <- exp(LnKmax[,"upr"]) |> round(digits = 3)
    kmax05 <- exp(LnKmax[,"lwr"]) |> round(digits = 3)
    kh95 <- (100*kmax95/PGA) |> round(digits = 1)
    kh05 <- (100*kmax05/PGA) |> round(digits = 1)
    KIT <- data.table(Da,kmax05,kmax,kmax95,kh05,kh,kh95,R2,RMSE,kmax_Units=TargetUnits)
  } else{
    KIT <- data.table(Da,kmax,kh,R2,RMSE,kmax_Units=TargetUnits)
  }
  return(KIT[kh>0 & kmax>0 & kmax<Inf])
}
