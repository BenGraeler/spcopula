# adopted from http://freakonometrics.hypotheses.org/2435, 

lowerEmpBivJointDepFun <- function(u) {
  stopifnot(ncol(u) == 2)
  empFun <- function(x) sum((u[,1]<=x)&(u[,2]<=x))/sum(u[,1]<=x)
  function(x) sapply(x,empFun)
}

upperEmpBivJointDepFun <- function(u) {
  stopifnot(ncol(u) == 2)
  empFun <- function(x) sum((u[,1]>=x)&(u[,2]>=x))/sum(u[,1]>=x)
  function(x) sapply(x,empFun)
}

empBivJointDepFun <- function(u) {
  stopifnot(ncol(u) == 2)
  
  function(z) {
    res <- z
    res[z>0.5] <- upperEmpBivJointDepFun(u)(z[z>0.5])
    res[z<=0.5] <- lowerEmpBivJointDepFun(u)(z[z<=0.5])
    return(res)
  }
}

##

lowerBivJointDepFun <- function(copula) {
  stopifnot(copula@dimension == 2)
  function(z) pCopula(cbind(z,z),copula)/z
}

upperBivJointDepFun <- function(copula) {
  stopifnot(copula@dimension == 2)
  function(z) (1-2*z+pCopula(cbind(z,z),copula))/(1-z)
}

bivJointDepFun <- function(copula) {
  function(z) {
    res <- z
    res[z>0.5] <- upperBivJointDepFun(copula)(z[z>0.5])
    res[z<=0.5] <- lowerBivJointDepFun(copula)(z[z<=0.5])
    return(res)
  }
}