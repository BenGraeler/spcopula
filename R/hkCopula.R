## Hierarchical Kendall Copulas as defined in Brechmann, Eike Christian. 
## "Hierarchical Kendall copulas: Properties and inference." Canadian Journal 
## of Statistics 42.1 (2014): 78-108.

# Slots:
#   
# Name:     nestingCop   clusterCops   kenFuns   dimension   parameters  param.names param.lowbnd  param.upbnd     fullname
# Class:        copula          list      list     integer      numeric    character      numeric      numeric    character

## nestingCop copula
## clusterCop list of list (copula, ind)

# easy constructor

hkCopula <- function(nestingCop, clusterCops, kenFuns=NULL) {
  if (is.null(kenFuns)) {
    kenFuns <- lapply(clusterCops, function(copInd) getKendallDistr(copInd[[1]]))
  }
  
  new("hkCopula", 
      nestingCop = nestingCop, 
      clusterCops = clusterCops,
      kenFuns = kenFuns,
      dimension = as.integer(sum(sapply(clusterCops, function(x) x[[1]]@dimension))+nestingCop@dimension-length(clusterCops)),
      parameters = NA_real_,
      param.names =NA_character_,
      param.lowbnd = NA_real_,
      param.upbnd = NA_real_, 
      fullname = "Hierarchical Kendall Copula")
}

showHkCopula <- function(object) {
  cat(object@fullname, "\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Nesting copula:\n")
  show(object@nestingCop)
  cat("Cluster copulas:\n")
  for (i in 1:length(object@clusterCops)) {
    cmpCop <- object@clusterCops[[i]][[1]]
    cat("  ", describeCop(cmpCop, "very short"), 
        "of dimension", cmpCop@dimension,
        "for indices", object@clusterCops[[i]][[2]], "\n")
  }
}

setMethod("show", signature("hkCopula"), showHkCopula)

## density

dHkCop  <- function(u, copula, log=F, ...) {
  stopifnot(ncol(u) == copula@dimension)

  lik <- NULL
  kenVal <- NULL
  
  for (i in 1:length(copula@clusterCops)) {
    cop <- copula@clusterCops[[i]][[1]]
    ind <- copula@clusterCops[[i]][[2]]
    ken <- copula@kenFuns[[i]]
    
    lik <- cbind(lik, dCopula(u[, ind], cop, log=log))
    kenVal <- cbind(kenVal, ken(pCopula(u[, ind], cop)))
  }
  
  if (ncol(kenVal) < copula@nestingCop@dimension) {
    kenVal <- cbind(kenVal, u[, -sapply(copula@clusterCops, function(x) x[[2]])])
  }
  
  lik <- cbind(lik, dCopula(kenVal, copula@nestingCop, log=log))
  
  if (log)
    return(apply(lik, 1, sum))
  
  return(apply(lik, 1, prod))
}

setMethod("dCopula", signature("matrix", "hkCopula"), dHkCop)
setMethod("dCopula", signature("numeric", "hkCopula"), 
          function(u, copula, log, ...) dHkCop(matrix(u, ncol = copula@dimension), copula, log, ...))

rHkCop <- function(n, copula, ...) {
  smpl <- matrix(NA, n, copula@dimension)

  nestSmpl <- rCopula(n, copula@nestingCop)
  
  for (i in 1:length(copula@clusterCops)) {
    cop <- copula@clusterCops[[i]][[1]]
    ind <- copula@clusterCops[[i]][[2]]
    ken <- copula@kenFuns[[i]]
    invKen <- genInvKenFun(ken)
    
    smpl[,ind] <- rCopula_y(invKen(nestSmpl[,i]), cop)
  }
  
  if (ncol(nestSmpl) > length(copula@clusterCops)) {
    smpl[,-sapply(copula@clusterCops, function(x) x[[2]])] <- nestSmpl[, -c(1:length(copula@clusterCops))]
  }
  
  return(smpl)
}

setMethod(rCopula, signature = c("numeric","hkCopula"), rHkCop)

setMethod(pCopula, signature = c("numeric","hkCopula"), 
          function(u, copula, ...) stop("Please use an empirical representation (i.e. \"genEmpCop\" applied to a sample of this copula)."))

setMethod(pCopula, signature = c("matrix","hkCopula"), 
          function(u, copula, ...) stop("Please use an empirical representation (i.e. \"genEmpCop\" applied to a sample of this copula)."))


rCop_y <- function(y, copula, n=1, n.disc = 1e2) {
  stopifnot(copula@dimension == 2)
  n.y <- length(y)
  stopifnot(n.y == 1 | n == 1)

  smpl <- matrix(NA, n.y*n, 2)
  
  for (i in 1:n.y) { # i <- 1
    condVals <- seq(y[i], 1-(1-y[i])/n.disc, length.out = n.disc)
    uv <- qCopula_v(copula, rep(y[i], n.disc), condVals)
    uv <- rbind(uv, qCopula_u(copula, rep(y[i], n.disc), condVals))
    
    uv <- uv[order(uv[,1]),]
    
    dSeq <- cumsum(c(0, apply((uv[-nrow(uv),]-uv[-1,])^2, 1, function (x) sqrt(sum(x)))))
    probs <- dCopula(uv, copula)
    
    apFun <- approxfun(dSeq, probs, rule = 2)
    probCor <- integrate(apFun, 0, max(dSeq))$value
    
    rContour <- runif(n, 0, probCor)
    
    funAppConPoint <- function(rCont) {
      invCDFContour <- function(x) {
        abs(integrate(apFun, 0, x)$value - rCont)
      } 
      
      lContour <- optimise(invCDFContour, c(0, max(dSeq)))$minimum
      
      dSeqInt <- findInterval(lContour, dSeq)
      
      lSeq <- sqrt(sum((uv[dSeqInt,]-uv[dSeqInt+1,])^2))
      
      uv[dSeqInt,] + (lContour - dSeq[dSeqInt])/lSeq * (uv[dSeqInt+1,]-uv[dSeqInt,])
    }
    
    if (n == 1) {
      appConPoint <- funAppConPoint(rContour)
      
      if (appConPoint[1] > appConPoint[2]) {
        smpl[i,] <- qCopula_u(copula, y[i], appConPoint[1])
      } else {
        smpl[i,] <- qCopula_v(copula, y[i], appConPoint[2])
      }
    } else {
      appConPoint <- t(sapply(rContour, funAppConPoint))
      
      boolLower <- appConPoint[,1] > appConPoint[,2]
      smpl[boolLower,] <- qCopula_u(copula, rep(y, sum(boolLower)), appConPoint[boolLower, 1])
      smpl[!boolLower,] <- qCopula_v(copula, rep(y, sum(!boolLower)), appConPoint[!boolLower, 2])
    }
    
    # plot(uv, type="l", xlim=c(uv[dSeqInt+c(0,1)]+c(-1,1)/1000), asp=1)
    # points(uv[dSeqInt+c(0,1),], col=c("red", "purple"))
    # points(matrix(appConPoint, nrow = 1), col="green")
    # points(matrix(smpl, nrow = 1), col="green", pch=2)
  }

  return(smpl)  
}

setGeneric("rCopula_y", function(y, copula, n=1, n.disc=1e2) NULL)
setMethod("rCopula_y", signature("numeric", "copula"), rCop_y)

# ## attic
# hkCop <- new("hkCopula", nestingCop=normalCopula(0.6),
#              clusterCops=list(list(cop=frankCopula(3), ind=c(1,2))),
#              kenFuns = list(getKendallDistr(frankCopula(3))),
#              dimension = 3L,
#              parameters = NA_real_,
#              param.names ="",
#              param.lowbnd = NA_real_,
#              param.upbnd = NA_real_,
#              fullname = "Hierarchical Kendall Copula")
# 
# hkCop4D <- new("hkCopula", nestingCop=normalCopula(0.6),
#                clusterCops=list(list(cop=frankCopula(3), ind=c(1,2)),
#                                 list(cop=gumbelCopula(5), ind=c(3,4))),
#                kenFuns = list(getKendallDistr(frankCopula(3)), 
#                               getKendallDistr(gumbelCopula(5))),
#                dimension = 4L,
#                parameters = c(0),
#                param.names ="",
#                param.lowbnd = c(0),
#                param.upbnd = 0, 
#                fullname = "Hierarchical Kendall Copula")
# 
# rHkCop(10, hkCop)
# 
# smplRHkCop3D <- rCopula(100, hkCop)
# 
# library(rgl)
# plot3d(smplRHkCop3D)
# kenHKcop <- genEmpKenFun(hkCop, sample = smplRHkCop3D)
# 
# curve(kenHKcop)
# 
# showMethods("pCopula")
# 
# plot(rCopula_y(0.9, gumbelCopula(5), 100))
# plot(rCopula_y(0.9, tawn3pCopula(c(0.75,.25,5)), 100))
# plot(rCopula_y(0.9, normalCopula(0.4), 100), asp=1)
# points(rCopula_y(0.9, normalCopula(0.8), 100), asp=1, pch=2)
# points(rCopula_y(0.9, normalCopula(-0.8), 100), asp=1, pch=3)
# 
# points(rCopula_y(0.4, normalCopula(-0.3), 100), asp=1, pch=4)
# abline(1.9,-1, col="red")
# 
# contour(normalCopula(-0.3), pCopula, asp=1)
# 
# sum(dHkCop(rHkCop(10, hkCop), hkCop))/10
# sum(dHkCop(rHkCop(10, hkCop4D), hkCop4D))/10
# 
# sum(dHkCop(matrix(runif(3000), 1000), hkCop))/1000
# sum(dHkCop(matrix(runif(4000), 1000), hkCop4D))/1000
# 
# par(mfrow=c(2,1))
# hist(dHkCop(matrix(runif(4*1e5),ncol = 4), hkCop4D), n=4000, xlim=c(0,10))
# hist(dHkCop(matrix(runif(3*1e5),ncol = 3), hkCop), n=400, xlim=c(0,10))
# 
# sum(dHkCop(matrix(runif(4*1e5),ncol = 4), hkCop4D))/1e5
# sum(dHkCop(matrix(runif(3*1e5),ncol = 3), hkCop))/1e5