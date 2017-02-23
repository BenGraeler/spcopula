############################################
## the spatial vine copula with covariate ##
############################################

validStCoVarVineCopula <- function(object) {
  return(is.function(object@coVarCop) & validStCopula(object@stCop) & validObject(object@topCop))
}

setClass("stCoVarVineCopula", representation("copula", coVarCop="function", stCop="stCopula", topCop="copula"),
         validity = validStCoVarVineCopula, contains=list("copula"))

## constructor ##
#################

stCoVarVineCopula <- function(coVarCop, stCop, topCop) {
  stopifnot(is(stCop,"stCopula"))
  
  new("stCoVarVineCopula", dimension = as.integer(topCop@dimension+1),
      parameters=numeric(), param.names = character(), param.lowbnd = numeric(), 
      param.upbnd = numeric(), 
      fullname = paste("Spatio-temporal covariate vine copula family with 1 spatio-temporal tree."),
      coVarCop=coVarCop, stCop=stCop, topCop=topCop)
}

## show ##
##########

showStCoVarVineCopula <- function(object) {
  dim <- object@dimension
  cat(object@fullname, "\n")
  cat("Dimension: ", dim, "\n")
}

setMethod("show", signature("stCoVarVineCopula"), showStCoVarVineCopula)

## density ##
#############

dstCoVarVine <- function(u, coVarCop, stCop, topCop, log, h, stInd) {
  stopifnot(nrow(u)==nrow(stInd))
  stopifnot(ncol(u)==2)
  stopifnot(dim(h)==3)

  l0 <- rep(0,nrow(u)) # level 0 spatio-temporal density
  dimDists <- dim(h)
  
  nrU <- nrow(u)
  
  u1 <- matrix(NA, nrU, ncol(u)-1)
  for(i in 2:dimDists[2]) { # i <- 1
    l0 <- l0 + dCopula(u[,c(1,i+1)], stCop, h=matrix(h[,i,], ncol=2), log=T)
    u1[,i] <- dduCopula(u[,c(1,i+1)], stCop, h=matrix(h[,i,], ncol=2))
  }
  
  uCoVar <- numeric(nrU)
  for (i in 1:nrU) {
    uCoVar[i] <- dduCopula(u[i,1:2], coVarCop(stInd[i,]))
  }
  
  if(!is.null(topCop))
    l1 <- dCopula(cbind(uCoVar, u1), topCop, log=T)
  else 
    l1 <- 0
  
  if(log)
    return(l0+l1)
  else(exp(l0+l1))
}

setMethod("dCopula", signature=signature("matrix","stCoVarVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dstCoVarVine(u, copula@coVarCop, copula@stCop, copula@topCop, log=log, ...)
            else
              dstCoVarVine(u, copula@coVarCop, copula@stCop, NULL, log=log, ...)
          })

setMethod("dCopula",signature=signature("numeric","stCoVarVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dstCoVarVine(matrix(u,ncol=copula@dimension), copula@coVarCop, copula@stCop, copula@topCop, log=log, ...)
            else
              dstCoVarVine(matrix(u,ncol=copula@dimension), copula@coVarCop, copula@stCop, NULL, log=log, ...)
          })

setMethod("dCopula",signature=signature("data.frame","stCoVarVineCopula"),
          function(u, copula, log, ...) {
            if("topCop" %in% slotNames(copula))
              dstCoVarVine(as.matrix(u), copula@stCop, copula@coVarCop, copula@stCop, copula@topCop, log=log, ...)
            else
              dstCoVarVine(as.matrix(u), copula@coVarCop, copula@stCop, NULL, log=log, ...)
          })
# 
# # fitting the spatial vine for a given list of spatial copulas
# fitStVine <- function(copula, data, method, estimate.variance=F) {
#   stopifnot(class(data)=="stNeighbourhood")
#   stopifnot(copula@dimension == ncol(data@data))
#   
#   u0 <- as.matrix(data@data) # previous level's (conditional) data
#   h0 <- data@distances # previous level's distances
#   l0 <- rep(0,nrow(u0)) # spatial density
#   u1 <- NULL # current level of conditional data
#   cat("[Margin ")
#   for(i in 1:dim(h0)[2]) { # i <- 1
#     l0 <- l0 + dCopula(u0[,c(1,i+1)], copula@stCop, h=h0[,i,], log=T)
#     cat(i,", ", sep="")
#     u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], copula@stCop, h=h0[,i,]))
#   }
#   u0 <- u1
#   cat("]\n")
#   
#   cat("[Estimating a",ncol(u0),"dimensional copula at the top.]\n")
#   vineCopFit <- fitCopula(copula@topCop, u0, method) 
#     
#   stVineCop <- stVineCopula(copula@stCop, vineCopFit@copula)
#   loglik <- vineCopFit@loglik
#   
#   return(new("fitCopula", estimate = stVineCop@parameters, var.est = matrix(NA), 
#              method = method, 
#              loglik = sum(l0)+loglik,
#              fitting.stats=list(convergence = as.integer(NA)),
#              nsample = nrow(data@data), copula=stVineCop))
# }

# setMethod("fitCopula",signature=signature("stVineCopula"),fitStVine)

# conditional spatio-temporal covariate vine
condStCoVarVine <- function (condVar, dists, stCVVC, stInd, n = 1000) {
  stopifnot(is.array(dists))
  
  coVarCop <- stCVVC@coVarCop(stInd)
  stBiCop <- stCVVC@stCop
  topCop <- stCVVC@topCop
  
  # add some points in the tails
  rat <- 50:1%x%c(1e-6,1e-5,1e-4,1e-3)
  xVals <- unique(sort(c(rat, 1 - rat, 1:(n - 1)/n)))
  nx <- length(xVals)
  nbs <- dim(dists)[2]
  
  repCondVar <- matrix(condVar, ncol = length(condVar), nrow = nx, byrow = T)
  dCoVar <- dCopula(cbind(xVals, repCondVar[,1]), coVarCop)
  condCoVar <- dduCopula(cbind(xVals, repCondVar[,1]), coVarCop)
  
  u0 <- cbind(xVals,repCondVar[,-1]) # previous level's (conditional) data
  condData <- matrix(NA, nx, nbs) # current level of conditional data
  dStBiC <- 1 # current likelihood
  
  for(i in 1:nbs) { # i <- 1
    condData[,i] <- dduCopula(u0[,c(1,i+1)], stBiCop, h=matrix(dists[,i,],1))
    dStBiC <- dStBiC*dCopula(u0[,c(1,i+1)], stBiCop, h=matrix(dists[,i,],1))
  }
  
  dTopVi <- dCopula(cbind(condCoVar,condData), topCop)
  
  density <- dCoVar*dStBiC*dTopVi
  # the 1-e6 corners linearily to [0,1], but ensure non-negative
  density <- c(max(0,2*density[1]-density[2]),
               density, max(0,2*density[nx]-density[nx-1]))
  linAppr <- approxfun(c(0, xVals, 1), density)
  
  # sum up the denstiy to rescale
  int <- sum(diff(c(0,xVals,1))*(0.5*diff(density)+density[-(nx+2)]))
  condVineFun <- function(u) linAppr(u)/int
  attr(condVineFun,"xVals") <- c(0,xVals,1)
  return(condVineFun)
}

# condition the covariate on the central location
condCovariate <- function(stNeigh, coVarCop) {
  stopifnot(length(stNeigh@coVar) == 1)
  
  nrNeigh <- nrow(stNeigh@index)
  ncData <- ncol(stNeigh@data)
  
  vddu <- numeric(nrNeigh)
  uv <- as.matrix(stNeigh@data[,c(1,ncData)])
  stInd <- stNeigh@index[,1,]
  
  for (i in 1:nrNeigh) {
    vddu[i] <- dduCopula(uv[i,], coVarCop(stInd[i,]))
  }
  
  return(vddu)
}


## interpolation ##
###################
# 
# stCopPredict.expectation <- function(predNeigh, stVine, margin, ..., stop.on.error=F) {
#   stopifnot(is.function(margin$q))
#   dists <- predNeigh@distances
#   
#   predMean <- NULL
#   for(i in 1:nrow(predNeigh@data)) { # i <-1
#     cat("[Predicting location ",i,".]\n", sep="")
#     condSecVine <- condStVine(as.numeric(predNeigh@data[i,]), dists[i,], stVine)
#     
#     condExp <-  function(x) {
#       margin$q(x)*condSecVine(x)
#     }
#     
#     ePred <- integrate(condExp,0,1,subdivisions=10000L,stop.on.error=stop.on.error, ...)
#     if(ePred$abs.error > 0.01)
#             warning("Numerical integration in predExpectation performed at a level of absolute error of only ",
#                     ePred$abs.error, " for location ",i,".")
#     predMean <- c(predMean, ePred$value)
#   }
#   if ("data" %in% slotNames(predNeigh@locations)) {
#     res <- predNeigh@locations
#     res@data[["expect"]] <- predMean
#     return(res)
#   } else {
#     predMean <- data.frame(predMean)
#     colnames(predMean) <- "expect"
#     return(addAttrToGeom(predNeigh@locations, predMean, match.ID=FALSE))
#   }
# }
# 
# stCopPredict.quantile <- function(predNeigh, stVine, margin, p=0.5) {
#   stopifnot(is.function(margin$q))
#   dists <- predNeigh@distances
#   
#   predQuantile <- NULL
#   for(i in 1:nrow(predNeigh@data)) { # i <-1
#     cat("[Predicting location ",i,".]\n", sep="")
#     condSecVine <- condStVine(as.numeric(predNeigh@data[i,]), dists[i,,,drop=F], stVine)
#     
#     xVals <- attr(condSecVine,"xVals")
#     density <- condSecVine(xVals)
#     nx <- length(xVals)
#     int <- cumsum(c(0,diff(xVals)*(0.5*diff(density)+density[-nx])))
#     lower <- max(which(int <= p))
#     m <- (density[lower+1]-density[lower])/(xVals[lower+1]-xVals[lower])
#     b <- density[lower]
#     xRes <- -b/m+sign(m)*sqrt(b^2/m^2+2*(p-int[lower])/m)
#     
#     predQuantile <- c(predQuantile, margin$q(xVals[lower]+xRes))
#   }
#   
#   if ("data" %in% slotNames(predNeigh@locations)) {
#     res <- predNeigh@locations
#     res@data[[paste("quantile.",p,sep="")]] <- predQuantile
#     return(res)
#   } else {
#     predQuantile <- data.frame(predQuantile)
#     colnames(predQuantile) <- paste("quantile.",p,sep="")
#     return(addAttrToGeom(predNeigh@locations, predQuantile, match.ID=FALSE))
#   }
# }
# 
# stCopPredict <- function(predNeigh, stVine, margin, method="quantile", p=0.5, ...) {
#   switch(method,
#          quantile=stCopPredict.quantile(predNeigh, stVine, margin, p),
#          expectation=stCopPredict.expectation(predNeigh, stVine, margin, ...))
# }
# 
# dropStTree <- function(neigh, stCop) {
#     stopifnot(class(neigh)=="stNeighbourhood")
#     
#     u0 <- as.matrix(neigh@data) # previous level's (conditional) data
#     h0 <- neigh@distances # previous level's distances
#     u1 <- NULL # current level of conditional data
#     h1s <- NULL # upcoming distances
#     h1t <- NULL # upcoming distances
#     cat("[Margin ")
#     for(i in 1:dim(h0)[2]) { # i <- 1
#       cat(i,", ")
#       u1 <- cbind(u1, dduCopula(u0[,c(1,i+1)], stCop, h=h0[,i,]))
#       if (i < ncol(neigh@distances)) {
#         h1s <- cbind(h1s, apply(neigh@index[, c(1, i + 1),1], 1, 
#                                 function(x) spDists(neigh@locations@sp[x, ])[1, 2]))
#         h1t <- cbind(h1t, apply(neigh@index[, c(1, i + 1),2], 1, 
#                                 function(x) diff(x)))
#       }
#     }
#     h1 <- array(dim=c(dim(h1s),2))
#     h1[,,1] <- h1s
#     h1[,,2] <- h1t
#     
#     varSplit <- strsplit(neigh@var, "|", fixed = TRUE)[[1]]
#     cond <- suppressWarnings(as.numeric(varSplit[length(varSplit)]))
#     if (is.na(cond)) 
#       cond <- paste(neigh@var, "|0", sep = "")
#     else cond <- paste(neigh@var, cond + 1, sep = "")
#     return(stNeighbourhood(data=u1, distances=h1, STxDF=neigh@locations, 
#                            ST=neigh@dataLocs, index=neigh@index[, -1,], 
#                            prediction=neigh@prediction, var=cond))
# }