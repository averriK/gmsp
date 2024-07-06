#' Title
#'
#' @param ATL data.table
#' @param TargetUnits character
#' @param kh_values vector
#'
#' @return
#' @export
#'
#' @examples



get_Newmark <- function(ATL,TargetUnits="mm",kh_values=c(0.01,0.02,0.05, 0.10, 0.15, 0.20, 0.25, 0.30,0.35, 0.40,0.45, 0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.90,1.00,1.25,1.50,1.75,2.00,2.50,3.00)){
  on.exit(expr={rm(list = ls())}, add = TRUE)
  . <- NULL
  g <- .getG(TargetUnits) #GMSP$g
  
  
  
  # Newmark Displacements  -------------------------------------------------------------
  
  NDL <- rbindlist(
    sapply(kh_values, function(kh) {
      ATL[, .(
        ID = sprintf("DN%02d", round(kh * 100)), 
        value = .get_ND(AT = s, t = t, kh = kh, FULL = FALSE)), 
        by = .(RecordSN, DIR, OCID)]}, simplify = FALSE))
  
  
  
  # Build Intensity Table
  NDW <- dcast(NDL, RecordSN + OCID + DIR ~ ID, value.var = "value")
  return(list(NDL=NDL,NDW=NDW))
}



.get_ND <- function(AT,t=NULL,dt=NULL,kh,TOL=1e-3, g=9806.650,FULL=TRUE) {
  on.exit(expr={rm(list = ls())}, add = TRUE)
  if(is.null(AT)) {return(NULL)}
  PGA <- max(abs(AT))
  if(PGA==0) {return(NULL)}
  ky <- kh*PGA/g
  ag <- AT/g
  NP <- length(ag)
  if(is.null(t) & !is.null(dt)){
    t <- seq(0,(NP-1)*dt,by=dt)
  }
  if(is.null(dt) & !is.null(t)) {
    dt <- mean(diff(t))
  } 
  
  a <- double(NP+1)
  v <- double(NP+1)
  u <- double(NP+1)
  for (i in 1:NP){
    if (v[i]<TOL){
      if (abs(ag[i])>ky){
        n <- ag[i]/abs(ag[i])
      }
      else {
        n <- ag[i]/ky
      }
    }
    else {
      n <- 1
    }
    a[i+1] <- (ag[i]-n*ky)*g
    # V(i+1)=V(i)+(1/2)*Dt*A(i)+1/2*Dt*A(i+1);
    v[i+1] <- v[i]+1/2*dt*(a[i+1]+a[i]) #ok
    if (v[i+1]<TOL){
      v[i+1] <- 0
      a[i+1] <- 0
    }
    # U(i+1)=U(i)+Dt*V(i)+(0.5-B)*(Dt^2)*A(i)+B*(Dt^2)*A(i+1);
    u[i+1] <- u[i]+dt*v[i]+(1/3)*(dt^2)*a[i]+1/6*(dt^2)*a[i+1] # Linear
    # u[i+1] <- u[i]+1/2*dt*(v[i+1]+v[i])# Constant
  }
  
  # return(Umax)
  if(FULL){
    return(u[1:NP])
  }
  else {
    Umax <- tail(u,1)
    return(Umax)
  }
}

