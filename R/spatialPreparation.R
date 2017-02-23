########################################################
##                                                    ##
## functions based on sp preparing the use of copulas ##
##                                                    ##
########################################################

## spatial neighbourhood constructor
####################################

neighbourhood <- function(data, distances, index, var, coVar=character(), prediction=FALSE) {
  sizeN <- ncol(distances)+1
  data <- as.data.frame(data)
  if (anyDuplicated(rownames(data))>0)
    rownames <- 1:length(rownames)
  
  new("neighbourhood", data=data, distances=distances, index=index,
      var=var, coVar=coVar, prediction=prediction)
}

## show
showNeighbourhood <- function(object){
  cat("A set of neighbourhoods consisting of", ncol(object@distances)+1, "locations each \n")
  if (length(object@var)>0) {
    cat("with", nrow(object@data), "rows of observations for:\n")
    cat(object@var, "\n")
  } else {
    cat("without data \n")
  }
  if(length(object@coVar)>0)
    cat("with covariate", object@coVar, "\n")
}

setMethod(show, signature("neighbourhood"), showNeighbourhood)

## names (from sp)
setMethod(names, signature("neighbourhood"), function(x) c(x@var,x@coVar))

selectFromNeighbourhood <- function(x, i) {
  new("neighbourhood", data=x@data[i,,drop=F], 
      distances=x@distances[i,,drop=F], index=x@index[i,,drop=F], 
      var=x@var, coVar=x@coVar, prediction=x@prediction)
}

setMethod("[", signature("neighbourhood","numeric"), selectFromNeighbourhood) 

## calculate neighbourhood from SpatialPointsDataFrame
getNeighbours <- function (dataLocs, predLocs, size = 5, 
                           var = names(dataLocs)[1], coVar=character(),
                           prediction = FALSE, min.dist = 0.01) {
  stopifnot((!prediction && missing(predLocs)) || (prediction && !missing(predLocs)))
  stopifnot(min.dist > 0 || prediction)
  
  if (missing(predLocs) && !prediction) 
    predLocs = dataLocs
  
  stopifnot(is(predLocs, "Spatial"))
  
  if ("data" %in% slotNames(dataLocs)) {
    if (any(is.na(match(var, names(dataLocs))))) 
      stop("The variables is not part of the data.")
  }
  
  nLocs <- length(predLocs)
  size <- min(size, length(dataLocs) + prediction)
  
  allLocs <- matrix(0, nLocs, size)
  allDists <- matrix(0, nLocs, size - 1)
  result = .C("getNeighbours_2d", 
     allLocs = as.double(allLocs),
     allDists = as.double(allDists),
     as.double(coordinates(dataLocs)),
     as.double(coordinates(predLocs)),
     as.integer(length(dataLocs)),
     as.integer(length(predLocs)),
     as.double(min.dist),
     as.integer(size), 
     as.integer(prediction),
     PACKAGE="spcopula")
  
  if ("data" %in% slotNames(dataLocs)) {
    if (!prediction) {
      allData <- matrix(dataLocs[result$allLocs, var, drop = F]@data[[1]],
                        nLocs, size)
    } else {
      allData <- matrix(c(rep(NA,nLocs),
                          dataLocs[result$allLocs[(nLocs+1):(nLocs*size)], var, drop = F]@data[[1]]),
                        nLocs, size)
    }
    colnames(allData) <- paste(paste("N", rep(0:(size - 1), each = length(var)), sep = ""),
                               rep(var, size), sep = ".")
  } else {
    allData <- as.data.frame(matrix(NA, nLocs, size + length(coVar)))
    var <- character()
  }
  
  return(neighbourhood(data=allData, distances=matrix(result$allDists, nLocs, size - 1), 
                       index=matrix(result$allLocs, nLocs, size), var=var, coVar=coVar,
                       prediction=prediction))
}

#############
## BINNING ##
#############

# calculates lag indicies for a Spatial object and stores the respective separating distances
# 
# boundaries  -> are the right-side limits in the dimenssion as provided by spDists
# data --------> a spatial object that can be handled by spDists()        
calcSpLagInd <- function(data, boundaries) {
  lags <- vector("list",length(boundaries))
  
  dists <- spDists(data)
  nlocs <- length(data)
  
  for (i in 1:(nlocs-1)) {
    for (j in (i+1):nlocs) {
      d <- dists[i,j]
      for ( k in 1:length(boundaries)) {
        if (d < boundaries[k]) {
          lags[[k]] <- rbind(lags[[k]],c(i,j,d))
          break()
        }
      }
    }
  }
  return(lags)
}

# the generic calcBins, calculates bins for spatial and spatio-temporal data
setGeneric("calcBins", function(data, var, nbins=15, boundaries=NA, cutoff=NA,
                                ..., cor.method="fasttau", plot=TRUE) {
                         standardGeneric("calcBins") 
                         })

## calculating the spatial bins
################################

calcSpBins <- function(data, var, nbins=15, boundaries=NA, cutoff=NA, ...,
                       cor.method="fasttau", plot=TRUE) {

  if(is.na(cutoff)) {
    cutoff <- spDists(coordinates(t(data@bbox)))[1,2]/3
  }
  if(any(is.na(boundaries))) {
    boundaries <- ((1:nbins) * cutoff/nbins)
  }
    
  nbins <- length(boundaries)-1
  
  lags <- calcSpLagInd(data, boundaries)
    
  mDists <- sapply(lags, function(x) mean(x[,3]))
  np <- sapply(lags, function(x) length(x[,3]))
  lagData <- lapply(lags, function(x) as.matrix((cbind(data[x[,1],var,drop=FALSE]@data, data[x[,2],var,drop=FALSE]@data))))
  
  if(cor.method == "fasttau")
    lagCor <- sapply(lagData, function(x) TauMatrix(x)[1,2])
  if(cor.method %in% c("kendall","spearman","pearson"))
    lagCor <- sapply(lagData, function(x) cor(x,method=cor.method)[1,2])
  if(cor.method == "normVariogram")  
    lagCor <- sapply(lagData, function(x) 1-cor(x,method="pearson")[1,2])
  if(cor.method == "variogram")  
    lagCor <- sapply(lagData, function(x) 0.5*mean((x[,1]-x[,2])^2,na.rm=T))
    
  if(plot) { 
    plot(mDists, lagCor, xlab="distance",ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(lagCor)),max(lagCor)), xlim=c(0,max(mDists)))
    abline(h=c(-min(lagCor),0,min(lagCor)),col="grey")
  }
  
#   res <- list(np=np, meanDists = mDists, lagCor=lagCor, lagData=lagData, lags=lags)
  res <- list(np=np, meanDists = mDists, lagCor=lagCor, lags=lags)
  attr(res,"cor.method") <- cor.method
  attr(res,"variable") <- var
  return(res)
}

setMethod(calcBins, signature("Spatial"), calcSpBins)

# calc bins from a (conditional) neighbourhood

calcNeighBins <- function(data, var=data@var, nbins=9, boundaries=NA, 
                          cutoff=NA, cor.method="kendall", plot=TRUE) {
  dists <- data@distances
  
  corFun <- switch(cor.method,
                   fasttau=function(x) TauMatrix(x)[1,2],
                   function(x) cor(x,method=cor.method)[1,2])
  
  if (any(is.na(boundaries))) 
    boundaries <- quantile(as.vector(dists), probs=c(1:nbins/nbins))
  if(!is.na(cutoff)) {
    boundaries <- boundaries[boundaries < cutoff]
    boundaries <- unique(c(0,boundaries,cutoff))
  } else {
    boundaries <- unique(c(0,boundaries))
  }
  
  nbins <- length(boundaries)-1
  
  np <- numeric(nbins)
  moa <- numeric(nbins)
  meanDists <- numeric(nbins)

  data <- as.matrix(data@data)
  
  lagData <- list()
  
  for (i in 1:nbins) {
    bools <- (dists <= boundaries[i+1] & dists > boundaries[i])
    
    pairs <- NULL
    for(col in 1:(dim(bools)[2])) {
      pairs <- rbind(pairs, data[bools[,col],c(1,1+col)])
    }
    
    lagData[[i]] <- pairs
    moa[i] <- corFun(pairs)
    meanDists[i] <- mean(dists[bools])
    np[i] <- sum(bools)
  }
  
  if(plot) { 
    plot(meanDists, moa, xlab="distance", ylab=paste("correlation [",cor.method,"]",sep=""), 
         ylim=1.05*c(-abs(min(moa, na.rm=T)),max(moa, na.rm=T)), xlim=c(0,max(meanDists,na.rm=T)))
    abline(h=c(-min(moa),0,min(moa)),col="grey")
  }
  
  res <- list(np=np, meanDists = meanDists, lagCor=moa, lagData=lagData)
  attr(res,"cor.method") <- switch(cor.method, fasttau="kendall", cor.method)
  attr(res,"variable") <- var
  
  return(res)
}
  
setMethod(calcBins, signature="neighbourhood", calcNeighBins)