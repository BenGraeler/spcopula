## inverse partial derivatives 
# numerical standard function
invdduCopula <- function(u, copula, y, ..., tol=.Machine$double.eps^0.5) {
    if (length(u) != length(y)) 
        stop("Length of u and y differ!")
    message("Numerical evaluation of invddu takes place.")
    res <- NULL
    for (i in 1:length(u)) {
        res <- rbind(res, optimize( function(x) 
          (dduCopula(cbind(rep(u[i], length(x)), x),copula) - y[i])^2, 
            interval = c(0, 1), tol=tol)$minimum)
    }
    return(res)
}

setGeneric("invdduCopula")

invddvCopula <- function(v, copula, y, ..., tol=.Machine$double.eps^0.5) {
    if (length(v) != length(y)) 
        stop("Length of v and y differ!")
  message("Numerical evaluation of invddv takes place.")
    res <- NULL
    for (i in 1:length(v)) {
        res <- rbind(res, optimize(function(x) 
          (ddvCopula(cbind(x, rep(v[i], length(x))),copula) - y[i])^2, 
            interval = c(0, 1), tol=tol)$minimum)
    }
    return(res)
}

setGeneric("invddvCopula")

###################
## Normal Copula ##
###################

## partial derivative d/du
##########################

dduNorm <- function(u, copula){
  rho <- copula@parameters

  u1 <- qnorm(u[,1]) # u ~ N(0,1)
  u2 <- qnorm(u[,2]) # v ~ N(0,1)

  return(pnorm(u2,mean=rho*u1,sd=sqrt(1-rho^2)))
}

setMethod("dduCopula", signature("numeric","normalCopula"),
          function(u, copula, ...) {
            dduNorm(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","normalCopula"), dduNorm)

## inverse of the partial derivative d/du
#########################################

invdduNorm <- function(u, copula, y){
  rho <- copula@parameters
  return(pnorm(qnorm(y,mean=rho*qnorm(u),sd=sqrt(1-rho^2))))
}

setMethod("invdduCopula", signature("numeric","normalCopula","numeric"), invdduNorm)


## partial derivative d/dv
##########################

ddvNorm <- function(u, copula){
  rho <- copula@parameters

  u1 <- qnorm(u[,1])
  u2 <- qnorm(u[,2])

  return(pnorm(u1,mean=rho*u2,sd=sqrt(1-rho^2)))
}

setMethod("ddvCopula", signature("numeric","normalCopula"),
          function(u, copula, ...) {
            ddvNorm(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","normalCopula"), ddvNorm)

## inverse of the partial derivative d/dv
#########################################

invddvNorm <- function(v, copula, y){
  rho <- copula@parameters
  return(pnorm(qnorm(y,mean=rho*qnorm(v),sd=sqrt(1-rho^2))))
}

setMethod("invddvCopula", signature("numeric","normalCopula","numeric"), invddvNorm)

########################
## independent copula ##
########################

## Kendall's tau
################
setMethod("tau", signature("indepCopula"), function(copula, ...) return(0))

## Spearman's rho
#################
setMethod("rho", signature("indepCopula"), function(copula, ...) return(0))

## indepCopula as evCopula derivatives of A
###########################################
setMethod("dAdu", signature("indepCopula"), 
          function(copula, w) {
            data.frame(der1=rep(0, length(w)), der2=rep(0, length(w)))
          })

## partial derivative d/du
##########################
setMethod("dduCopula", signature("numeric","indepCopula"),
          function(u, copula, ...) {
            matrix(u,ncol=copula@dimension)[,2]
          })
setMethod("dduCopula", signature("matrix","indepCopula"), function(u, copula, ...) u[,2])

## inverse of the partial derivative d/du
#########################################
invdduIndep <- function(u, copula, y){
  return(y)
}

setMethod("invdduCopula", signature("numeric","indepCopula","numeric"), invdduIndep)

## partial derivative d/dv
##########################
setMethod("ddvCopula", signature("numeric","indepCopula"),
          function(u, copula, ...) {
            matrix(u,ncol=copula@dimension)[,1]
          })
setMethod("ddvCopula", signature("matrix","indepCopula"), function(u, copula, ...) u[,1])

## inverse of the partial derivative d/dv
#########################################
invddvIndep <- function(v, copula, y){
  return(y)
}

setMethod("invddvCopula", signature("numeric","indepCopula", "numeric"), invddvIndep)


####################
## Clayton Copula ##
####################

## partial derivative d/du
##########################
dduClayton <- function(u, copula){
  rho <- copula@parameters

  if (rho==0)
    return(u[,2])
  
  u1 <- u[,1]
  u2 <- u[,2]

  pmax(u1^(-rho)+u2^(-rho)-1,0)^((-1-rho)/rho)*u1^(-rho-1)
}

setMethod("dduCopula", signature("numeric","claytonCopula"),
          function(u, copula, ...) {
            dduClayton(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","claytonCopula"), dduClayton)

## inverse of the partial derivative d/du
#########################################

invdduClayton <- function(u, copula, y){
    rho <- copula@parameters[1]
    if (length(u)!=length(y)) 
        stop("Length of u and y differ!")
    return(((y^(rho/(-1-rho))-1)*u^(-rho)+1)^(-1/rho)) # by DL
}

setMethod("invdduCopula", signature("numeric","claytonCopula","numeric"), invdduClayton)

## partial derivative d/dv
##########################

ddvClayton <- function(u, copula){
  rho <- copula@parameters
  if (!is.matrix(u)) u <- matrix(u, ncol=2)

  if (rho==0)
    return(u[,1])
  
  u1 <- u[,1]
  u2 <- u[,2]

  pmax(u2^(-rho)+u1^(-rho)-1,0)^((-1-rho)/rho)*u2^(-rho-1)
}

setMethod("ddvCopula", signature("numeric","claytonCopula"),
          function(u, copula, ...) {
            ddvClayton(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","claytonCopula"), ddvClayton)


## inverse of the partial derivative d/dv
#########################################

invddvClayton <- function(v, copula, y){
    rho <- copula@parameters[1]
    if (length(v)!=length(y)) 
        stop("Length of v and y differ!")
    return(((y^(rho/(-1-rho))-1)*v^(-rho)+1)^(-1/rho))
}

setMethod("invddvCopula", signature("numeric", "claytonCopula", "numeric"), invddvClayton)


###################
## Gumbel Copula ##
###################

## partial derivative d/du
##########################

dduGumbel <- function(u, copula){
  rho <- copula@parameters

  u1 <- u[,1]
  u2 <- u[,2]

  pCopula(u,gumbelCopula(rho)) * ((-log(u1))^rho+(-log(u2))^rho)^(1/rho-1) * (-log(u1))^(rho-1)/u1
}

setMethod("dduCopula", signature("numeric","gumbelCopula"),
          function(u, copula, ...) {
            dduGumbel(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","gumbelCopula"), dduGumbel)


## partial derivative d/dv
##########################

ddvGumbel <- function(u, copula){
  rho <- copula@parameters

  u1 <- u[,1]
  u2 <- u[,2]

  pCopula(u,gumbelCopula(rho)) * ((-log(u2))^rho+(-log(u1))^rho)^(1/rho-1) * (-log(u2))^(rho-1)/u2
}

setMethod("ddvCopula", signature("numeric","gumbelCopula"),
          function(u, copula, ...) {
            ddvGumbel(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","gumbelCopula"), ddvGumbel)



##################
## Frank Copula ##
##################

## partial derivative d/du
##########################

## Wolfram alpha:
# -e^a/(e^(a u) - e^a) - ((e^a - 1) e^(a + a u))/((e^(a u) - e^a) (e^a - e^(a + a u) + e^(a u + a v) - e^(a + a v)))

dduFrank <- function(u, copula){
  rho <- copula@parameters

  v <- u[,2]
  u <- u[,1]
  E <- exp(1)
  eRho <- E^rho
  eRhoU <- E^(rho * u)
  eRhoV <- E^(rho * v)
  
  -eRho/(eRhoU - eRho) - ((eRho - 1) * eRho * eRhoU)/((eRhoU - eRho) * (eRho - eRho * eRhoU + eRhoU * eRhoV - eRho * eRhoV))
}

setMethod("dduCopula", signature("numeric","frankCopula"),
          function(u, copula, ...) {
            dduFrank(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","frankCopula"), dduFrank)

## partial derivative d/dv
##########################

ddvFrank <- function(u, copula){
  rho <- copula@parameters
  
  v <- u[,2]
  u <- u[,1]
  E <- exp(1)
  eRho <- E^rho
  eRhoU <- E^(rho * u)
  eRhoV <- E^(rho * v)
  
  -eRho/(eRhoV - eRho) - ((eRho - 1) * eRho * eRhoV)/((eRhoV - eRho) * (eRho - eRho * eRhoV + eRhoV * eRhoU - eRho * eRhoU))
}

setMethod("ddvCopula", signature("numeric","frankCopula"),
          function(u, copula, ...) {
            ddvFrank(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","frankCopula"), ddvFrank)

####################
## student Copula ##
####################

## partial derivative d/du
##########################

dduStudent <- function(u, copula){
  df <- copula@parameters[2]
  v <- qt(u,df=df)
  
  rho <- copula@parameters[1]
  
  return(pt(sqrt((df+1)/(df+v[,1]^2)) / sqrt(1 - rho^2) * (v[,2] - rho * v[,1]), df=df+1))
}

setMethod("dduCopula", signature("numeric","tCopula"),
          function(u, copula, ...) {
            dduStudent(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("dduCopula", signature("matrix","tCopula"), dduStudent)


## partial derivative d/dv
##########################

ddvStudent <- function(u, copula){
  df <- copula@parameters[2]
  v <- qt(u, df=df)
  
  rho <- copula@parameters[1]
  
  return(pt(sqrt((df+1)/(df+v[,2]^2)) / sqrt(1 - rho^2) * (v[,1] - rho * v[,2]), df=df+1))
}

setMethod("ddvCopula", signature("numeric","tCopula"),
          function(u, copula, ...) {
            ddvStudent(matrix(u,ncol=copula@dimension),copula)
          })
setMethod("ddvCopula", signature("matrix","tCopula"), ddvStudent)