#######################################
## tawn copula with all 3 parameters ##
#######################################

setClass("tawn3pCopula", representation(exprdist = "expression"),
         contains = "evCopula")

Atawn3p <- function(t, param = c(0.9302082, 1, 8.355008)) {
  alpha <- param[1]
  beta  <- param[2]
  theta <- param[3]
  (1-beta)*(t) + (1-alpha)*(1-t) + ((alpha*(1-t))^theta+(beta*t)^theta)^(1/theta)

}

ATawn <- function(copula, w) {
  Atawn3p(w, copula@parameters)
}

setMethod("A", signature("tawn3pCopula"), ATawn)

dAduTawn <- function(copula, w) {
  alpha <- copula@parameters[1]
  beta  <- copula@parameters[2]
  theta <- copula@parameters[3]

  # 1st derivative
  p1 <- (alpha*(alpha*(-(w-1)))^(theta-1)-beta*(beta*w)^(theta-1)) 
  p2 <- ((alpha*(-(w-1)))^theta+(beta*w)^theta)^(1/theta-1)
  
  # 2nd derivative
  p3 <- (alpha*(-(w-1)))^(theta-2)
  p4 <- (beta*w)^(theta-2)
  p5 <- ((alpha*(-(w-1)))^theta+(beta*w)^theta)^(1/theta-2)
  
  data.frame(der1=alpha-beta-p1*p2,
             der2=alpha^2*beta^2*(theta-1)*p3*p4*p5)
}

setMethod("dAdu", signature("tawn3pCopula"), dAduTawn)

tawn3pCopula <- function (param = c(0.5, 0.5, 2)) {
  # A(t) = (1-beta)*t + (1-alpha)*(1-t) + ((alpha*(1-t))^theta+(beta*t)^theta)^(1/theta)
  # C(u1,u2) = exp(log(u1*u2) * A(log(u2)/log(u1*u2)))
  #          = u1*u2 + exp(A(log(u2)/log(u1*u2)))

  cdf <- expression(exp(log(u1*u2)*((1-beta)*(log(u2)/log(u1*u2)) +
                                    (1-alpha)*(1-log(u2)/log(u1*u2)) +
                                    ((alpha*(1-log(u2)/log(u1*u2)))^theta+(beta*log(u2)/log(u1*u2))^theta)^(1/theta))))
  dCdU1 <- D(cdf, "u1")
  dCdU2 <- D(cdf, "u2")
  pdf <- D(dCdU1, "u2")
  
  new("tawn3pCopula", dimension = 2L, exprdist = c(cdf = cdf, pdf = pdf,
                                                   dCdU = dCdU1, dCdV = dCdU2),
      parameters = param, param.names = c("alpha", "beta", "theta"), 
      param.lowbnd = c(0,0,1), param.upbnd = c(1,1,Inf), 
      fullname = "Tawn copula family with three parameters; Extreme value copula")
}

dtawn3pCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  for (i in 1:dim) {
    assign(paste("u", i, sep=""), u[,i])
  }
  alpha <- copula@parameters[1]
  beta  <- copula@parameters[2]
  theta <- copula@parameters[3]
  
  val <- c(eval(copula@exprdist$pdf))
  ## improve log-case
  if(log) 
    return(log(val))
  else 
    val
}

setMethod("dCopula", signature(copula = "tawn3pCopula"), dtawn3pCopula)

ptawn3pCopula <- function(u, copula, ...) {
  dim <- copula@dimension
  for (i in 1:dim) {
    assign(paste("u", i, sep=""), u[,i])
  }
  alpha <- copula@parameters[1]
  beta <-  copula@parameters[2]
  theta <-copula@parameters[3]
  
  val <- c(eval(copula@exprdist$cdf))
}

setMethod("pCopula", signature(copula = "tawn3pCopula"),  ptawn3pCopula)

# partial derivatives

ddutawn3pCopula <- function(u, copula, ...) {
  dim <- copula@dimension
  for (i in 1:dim) {
    assign(paste("u", i, sep=""), u[,i])
  }
  
  alpha <- copula@parameters[1]
  beta  <- copula@parameters[2]
  theta <- copula@parameters[3]
  
  return(eval(copula@exprdist$dCdU))
}

setMethod("dduCopula", signature(copula = "tawn3pCopula"), ddutawn3pCopula)

ddvtawn3pCopula <- function(u, copula, ...) {
  dim <- copula@dimension
  for (i in 1:dim) {
    assign(paste("u", i, sep=""), u[,i])
  }
  
  alpha <- copula@parameters[1]
  beta  <- copula@parameters[2]
  theta <- copula@parameters[3]
  
  return(eval(copula@exprdist$dCdV))
}

setMethod("ddvCopula", signature(copula = "tawn3pCopula"), ddvtawn3pCopula)

# tawn3pCop <- tawn3pCopula()
# dduCopula(cbind(runif(10), runif(10)), tawn3pCop)
## fit

# fitTawn3pCop <- function(copula, data, method = c("mpl", "ml"), 
#                          start = copula@parameters,
#                          lower = copula@param.lowbnd,
#                          upper = copula@param.upbnd,
#                          optim.method = "L-BFGS-B", 
#                          optim.control = list(maxit = 1000), estimate.variance = FALSE, 
#                          hideWarnings = TRUE) {
#   
#   fitCopulaAny <- selectMethod(fitCopula, "copula")
#   fitCopulaAny(copula, data, method, start, lower, upper, 
#                optim.method, optim.control, estimate.variance,
#                hideWarnings)
# }
#             
# setMethod("fitCopula", signature("tawn3pCopula"), fitTawn3pCop)