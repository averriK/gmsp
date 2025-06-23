#' get_ND
#'
#' @param AT data.table
#' @param t numeric
#' @param dt numeric
#' @param kh numeric
#' @param TOL numeric
#' @param g numeric
#' @param full boolean
#'
#' @return list
#' @export
#'

get_ND <- function(AT, t = NULL, dt = NULL, kh, TOL = 1e-3, g = 9806.650, full = TRUE) {
    on.exit(expr = {
        rm(list = ls())
    }, add = TRUE)
    if (is.null(AT)) {
        return(NULL)
    }
    PGA <- max(abs(AT))
    if (PGA == 0) {
        return(NULL)
    }
    ky <- kh * PGA / g
    ag <- AT / g
    NP <- length(ag)
    if (is.null(t) & !is.null(dt)) {
        t <- seq(0, (NP - 1) * dt, by = dt)
    }
    if (is.null(dt) & !is.null(t)) {
        dt <- mean(diff(t))
    }

    a <- double(NP + 1)
    v <- double(NP + 1)
    u <- double(NP + 1)
    for (i in 1:NP) {
        if (v[i] < TOL) {
            if (abs(ag[i]) > ky) {
                n <- ag[i] / abs(ag[i])
            } else {
                n <- ag[i] / ky
            }
        } else {
            n <- 1
        }
        a[i + 1] <- (ag[i] - n * ky) * g
        # V(i+1)=V(i)+(1/2)*Dt*A(i)+1/2*Dt*A(i+1);
        v[i + 1] <- v[i] + 1 / 2 * dt * (a[i + 1] + a[i]) # ok
        if (v[i + 1] < TOL) {
            v[i + 1] <- 0
            a[i + 1] <- 0
        }
        # U(i+1)=U(i)+Dt*V(i)+(0.5-B)*(Dt^2)*A(i)+B*(Dt^2)*A(i+1);
        u[i + 1] <- u[i] + dt * v[i] + (1 / 3) * (dt^2) * a[i] + 1 / 6 * (dt^2) * a[i + 1] # Linear
        # u[i+1] <- u[i]+1/2*dt*(v[i+1]+v[i])# Constant
    }

    # return(Umax)
    if (full == TRUE) {
        return(u[1:NP])
    } else {
        Umax <- tail(u, 1)
        return(Umax)
    }
}
