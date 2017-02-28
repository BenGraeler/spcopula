genEmpKenFun <- function(copula, sample=NULL) {
  if(is.null(sample)) 
    sample <- rCopula(1e6,copula)
  if(missing(copula)) {
    # taken from package copula function "Kn"
    stopifnot((n <- nrow(sample)) >= 1, (d <- ncol(sample)) >= 1)
    ken <- vapply(seq_len(n), function(i) sum(colSums(t(sample) < sample[i, ]) == d)/(n + 1), NA_real_)
  } else {
    ken <- pCopula(sample, copula)
  }
  
  return(ecdf(ken))
}

## inverse kendall function
# Kendall return period:
# K_x(t)= \mu / (1-K_C(t))
# 
# solve K_x(t)=1/1000
# 
# KRP = \mu / (1-K_C(t))
# <=> (1-K_C(t)) = \mu / KRP
# <=> K_C(t) = 1 - \mu /KRP
# <=> t = K_C^{-1}(1 - \mu /KRP)

genInvKenFun <- function(kenFun, tol=.Machine$double.eps^.5) {
  invKenFun <- function(k){
    res <- NULL
    for(i in 1:length(k)) {
      res <- c(res, optimize(function(x) (kenFun(x)-k[i])^2, c(0,1), tol=tol)$minimum)
    }
    return(res)
  }
  return(invKenFun)
}

## return periods
kendallRP <- function(kendallFun, cl=c(.99,.999), mu=1, copula) {
  if(missing(kendallFun) & missing(copula)) 
      stop("Either the kendall distribution function or the copula must be provided. Note that the calculation of the kendall distribution function from the copula is pretty time consuming. Saving them separately might be advantageous.")
  if(missing(kendallFun)) kendallFun <- genEmpKenFun(copula)
  if(length(mu)>1 & length(cl) > 1) stop("Either the critial level (cl) or mu may be of length larger than 1!")
  return(mu/(1-kendallFun(cl)))
}   

criticalLevel <- function(kendallFun, KRP=c(100,1000), mu=1, copula) {
  if(missing(kendallFun) & missing(copula)) 
      stop("Either the kendall distribution function or the copula must be provided. Note that the calculation of the kendall distribution function from the copula is pretty time consuming. Saving them separately might be advantageous.")
  if(missing(kendallFun))
      kendallFun <- genEmpKenFun(copula)
  if(length(mu)>1 & length(KRP) > 1) 
      stop("Either the kendall return period or mu may be of length larger than 1!")
  invKenFun <- genInvKenFun(kendallFun)
  return(invKenFun(1-mu/KRP))
}

## next: calculating critical layer, sampling from the layer, selecting "typical" points
# calculate critical layer (ONLY 2D by now)
criticalPair <- function(copula, cl, u, ind, tol=sqrt(.Machine$double.eps)) {
  
  optimFun <- function(x, u, ind) {
    pair <- cbind(x,x)
    pair[,ind] <- u
    return(abs(pCopula(pair,copula)-cl))
  }
  
  sapply(u, function(uRow) {
              upper <- cl+(1-uRow)
              if (upper<cl | uRow < cl) 
                return(NA)
              if (upper == cl) 
                return(cl)
              optimize(function(x) optimFun(x, uRow, ind),
                       interval=c(cl,upper), tol=tol)$minimum
            })
}


# calculate critical layer (ONLY 3D by now)
criticalTriple <- function(copula, cl, u, ind, tol=sqrt(.Machine$double.eps)) {
  if(!is.matrix(u)) u <- matrix(u,ncol=2)
    
  optimFun <- function(x, u, ind) {
    x <- matrix(rep(x,3),ncol=3)
    x[,ind[1]] <- u[1]
    x[,ind[2]] <- u[2]
    return(abs(pCopula(x,copula)-cl))
  }
  
  apply(u, 1, 
        function(uRow) {
          upper <- min(1,2+cl-sum(uRow)) # hyperplane in the hypercube
          if (upper < cl | any(uRow < cl)) 
            return(NA)
          if (upper == cl) 
            return(cl)
          optimize(function(x) optimFun(x, uRow, ind), 
                   interval=c(cl,upper),tol=tol)$minimum
        })
}

qCopula_u.def <- function(copula, p, u, tol=.Machine$double.eps^.5) { # sample=NULL
  copDim <- dim(copula)
  stopifnot(length(p) == length(u)) 
  
  if (copDim == 2) {
    res <- sapply(1:length(p), 
                  function(ind) {
                    if (u[ind] < p[ind]) 
                      return(NA)
                    if (u[ind] == 1)
                      return(p[ind])
                    optimise(function(v) abs(pCopula(cbind(u[ind], v), copula) - p[ind]),
                             c(p[ind], 1 + p[ind] - u[ind]), tol=tol)$minimum
                  })
  } else {
  res < NULL
    for(i in 1:length(p)) { # i <- 1
      if (u[i] < p[i]) {
        res <- rbind(res, rep(NA,dim-1))
      } else {
        opt <- optim(par=rep(p[i],dim-1), 
                     function(vw) abs(pCopula(c(u[i],vw), copula)-p[i]), 
                     lower=rep(p[i],dim-1), upper=rep(1,dim-1), method="L-BFGS-B")
        res <- rbind(res, opt$par)
      }
    }
  }
  
  return(cbind(u, res))
}

setGeneric("qCopula_u", function(copula, p, u, ...) standardGeneric("qCopula_u"))
setMethod("qCopula_u", signature("copula"), qCopula_u.def)


qCopula_v.def <- function(copula, p, v, tol=.Machine$double.eps^.5) {
  copDim <- dim(copula)
  if(length(p) != length(v)) 
    stop("Length of p and u differ!")
  
  if (copDim == 2) {
    res <- sapply(1:length(p), 
                  function(ind) {
                    if (v[ind] < p[ind]) 
                      return(NA)
                    if (v[ind] == 1)
                      return(p[ind])
                    optimise(function(u) abs(pCopula(cbind(u, v[ind]), copula) - p[ind]),
                             c(p[ind], 1 + p[ind] - v[ind]), tol=tol)$minimum
                  })
    res <- cbind(res, v)
  } else {
    res < NULL
    for(i in 1:length(p)) { # i <- 1
      if (v[i] < p[i]) {
        res <- rbind(res,rep(NA,dim-1))
      } else {
        opt <- optim(par=rep(p[i],dim-1), 
                     function(uw) abs(pCopula(c(uw[1], v[i], uw[2]), copula)-p[i]), 
                     lower=rep(p[i],dim-1), upper=rep(1,dim-1), method="L-BFGS-B")
        res <- rbind(res, opt$par)
      }
    }
    
    res <- cbind(res[,1], v, res[,2])
  }
  
  return(res)
}

setGeneric("qCopula_v", function(copula, p, v, ...) standardGeneric("qCopula_v"))
setMethod("qCopula_v", signature("copula"), qCopula_v.def)


## kendall distribution

# generic kendall function
kendall <- function(t, copula) {
  standardGeneric(t, copula)
}

# empirical default
getKendallDistr <- function(copula, sample=NULL) {
#   standardGeneric("getKendallDistr")
  if(is.null(sample))
    sample <- rCopula(1e6, copula)
  empCop <- empiricalCopula(sample, copula)
  ken <- pCopula(sample, empCop) # takes really long, any suggestions? Comparring a 1e6x3/1e6x2 matrix by 1e6 pairs/triplets values
  
  empKenFun <- function(tlevel) {
    res <- NULL
    for(t in tlevel) {
      res <- c(res,sum(ken<=t))
    }
    return(res/nrow(sample))
  }
  return(empKenFun)
}

setGeneric("getKendallDistr")

## 

kendallDistribution <- function(copula, t) {
  stop("There is no analytical expression implemented for this copula family. See 'getKendallDistr' for a numerical solution instead.")
}

setGeneric("kendallDistribution")

## Clayton
## kendall distribution/measure, taken from VineCopula:::obs.stat
kendall.Clayton <- function(copula, t){
  par = copula@parameters
  
  kt <- rep(NA,length(t))
  kt <- t + t * (1 - t^par)/par
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("claytonCopula"), kendall.Clayton) # for easy backwards compatibility
setMethod("kendall", signature("numeric", "claytonCopula"), function(t, copula) kendall.Clayton(copula, t))

setMethod("getKendallDistr", signature("claytonCopula"), 
          function(copula) return(function(t) kendall.Clayton(copula, t)))

## Gumbel
## kendall distribution/measure, taken from VineCopula:::obs.stat
kendall.Gumbel <- function(copula, t){
  par = copula@parameters
  
  kt <- rep(NA,length(t))
  kt <- t - t * log(t)/(par)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("gumbelCopula"), kendall.Gumbel) # for easy backwards compatibility
setMethod("kendall", signature("numeric","gumbelCopula"), function(t, copula) kendall.Gumbel(copula, t)) 

setMethod("getKendallDistr", signature("gumbelCopula"), 
          function(copula) return(function(t) kendall.Gumbel(copula, t)))

## Frank
## kendall distribution/measure, taken from VineCopula:::obs.stat
kendall.Frank <- function(copula, t){
  par = copula@parameters
  
  kt <- rep(NA,length(t))
  kt <- t + log((1 - exp(-par))/(1 - exp(-par * t))) * (1 - exp(-par * t))/(par * exp(-par * t))
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("frankCopula"), kendall.Frank) # for easy backwards compatibility
setMethod("kendall", signature("numeric", "frankCopula"),  function(t, copula) kendall.Frank)

setMethod("getKendallDistr", signature("frankCopula"), 
          function(copula) return(function(t) kendall.Frank(copula, t)))

## direct definition for Archimedean copulas

# BB1
## kendall distribution/measure
kendall.BB1 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + 1/(theta * delta) * (t^(-theta) - 1)/(t^(-1 - theta))
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB1Copula"), kendall.BB1)
setMethod("kendall", signature("numeric", "BB1Copula"), 
          function(t, copula) {
            stopifnot(copula@dimension <= 2) 
            kendall.BB1(copula, t)
          })

setMethod("getKendallDistr", signature("BB1Copula"), function(copula) return(function(t) kendall.BB1(copula, t)) )


# BB6
## kendall distribution/measure, taken from VineCopula:::obs.stat
kendall.BB6 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + log(-(1 - t)^theta + 1) * (1 - t - (1 - t)^(-theta) + (1 - t)^(-theta) * t)/(delta * theta)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB6Copula"), kendall.BB6) # for easy backwards compatibility
setMethod("kendall", signature("numeric", "BB6Copula"), 
          function(t, copula) {
            stopifnot(copula@dimension <= 2) 
            kendall.BB6(copula, t)
          })

setMethod("getKendallDistr", signature("BB6Copula"), 
          function(copula) return(function(t) kendall.BB6(copula, t)))

# BB7
## kendall distribution/measure, taken from VineCopula:::obs.stat
kendall.BB7 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + 1/(theta * delta) * ((1 - (1 - t)^theta)^(-delta) -  1)/
    ((1 - t)^(theta - 1) * (1 - (1 - t)^theta)^(-delta - 1))
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB7Copula"), kendall.BB7)
setMethod("kendall", signature("numeric", "BB7Copula"), 
          function(t, copula) {
            stopifnot(copula@dimension <= 2) 
            kendall.BB7(copula, t)
          })

setMethod("getKendallDistr", signature("BB7Copula"), 
          function(copula) return(function(t) kendall.BB7(copula, t)))

# BB8
## kendall distribution/measure, taken from VineCopula:::obs.stat
kendall.BB8 <- function(copula, t){
  theta = copula@parameters[1]
  delta = copula@parameters[2]
  
  kt <- rep(NA,length(t))
  kt <- t + log(((1 - t * delta)^theta - 1)/((1 - delta)^theta - 1)) * (1 - t * delta - (1 - t * delta)^(-theta) + (1 - t * delta)^(-theta) * t * delta)/ (theta * delta)
  kt[t==1] <- 1
  kt[t==0] <- 0
  return(kt)  
}

setMethod("kendallDistribution", signature("BB8Copula"), kendall.BB8) # for easy backwards compatibility
setMethod("kendall", signature("numeric", "BB8Copula"), 
          function(t, copula) {
            stopifnot(copula@dimension <= 2) 
            kendall.BB8(copula, t)
          })

setMethod("getKendallDistr", signature("BB8Copula"), 
          function(copula) return(function(t) kendall.BB8(copula, t)))

# BiJoe
kendall.Joe <- function(copula, t) kdJoe(t, copula)
  
  ## kendall distribution/measure, taken from VineCopula:::obs.stat
#   {
#   par = copula@parameters[1]
#   
#   kt <- rep(NA,length(t))
#   kt <- t - (log(1 - (1 - t)^par) * (1 - (1 - t))^par)/(par * (1 - t)^(par - 1))
#   kt[t==1] <- 1
#   kt[t==0] <- 0
#   return(kt)  
# }

setMethod("kendallDistribution", signature("joeBiCopula"), kendall.Joe)

setMethod("getKendallDistr", signature("joeBiCopula"), 
          function(copula) return(function(t) kendall.Joe(copula, t)))

