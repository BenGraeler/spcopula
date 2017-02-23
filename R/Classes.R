####
## an empirical survival copula representation

validEmpSurCopula <- function(object) {
  if(ncol(object@sample) != object@dimension)
    return("Dimension of the copula and the sample do not match.")
  else
    return(TRUE)
}

setClass("empSurCopula",
         representation = representation("copula", sample="matrix"),
         validity = validEmpSurCopula,
         contains = list("copula")
)

####
## the leaf copula

validLeafCopula <- function(object) {
  if (object@dimension != 2)
    return("The leaf copula only supports two dimensions.")
  
  if (any(is.na(object@parameters)))
    return("Parameter value is \"NA\".")
  else return (TRUE)
}

setClass("leafCopula",
         representation = representation("copula"),
         validity = validLeafCopula,
         contains = list("copula")
)

#####
## Hierarchical Kendall Copulas

validHkCopula <- function(object) {
  stopifnot(all(sapply(object@clusterCops, function(x) inherits(x[[1]], "copula"))))
  if (object@dimension != sum(sapply(object@clusterCops, function(x) x[[1]]@dimension))+object@nestingCop@dimension-length(object@clusterCops))
    return("The dimensions of the hierarchical Kendall copula do not match.")
  
  if(length(object@clusterCops) != length(object@kenFuns))
    return("Each cluster copula needs to have its Kendall function in 'kenFuns'.")
  
  else return (TRUE)
}

setClass("hkCopula",
         representation = representation("copula",
                                         nestingCop = "copula",
                                         clusterCops = "list",
                                         kenFuns = "list"),
         validity = validHkCopula,
         contains = list("copula")
)




## 
## the spatial copula
##
## realized as a distance dependent convex combination of biv copulas

# dimension = "numeric"     set to 2
# parameters = "numeric"    set of parameters
# param.names = "character" appropriate names
# param.lowbnd = "numeric"  appropriate lower bounds
# param.upbnd = "numeric"   appropriate upper bounds
# fullname = "character"    name printed with "show"
# components="list"         list of copulas 
# distances="numeric"       the linking distances
# unit="character"          measurement unit of distance
# depFun="function"         an optional dependence function; depFun(NULL)
#                             has to return either "spearman" or "kendall" 
#                             dependening on the moa used. Make sure depFun
#                             assings valid parameters to the copulas involved

validSpCopula <- function(object) {
  if (length(object@components) != length(object@distances)) 
    return("Length of components does not equal length of distances. \n Note: The last distance must give the range and it is automatically associated with the indepenence copula.")
  check.upper <- NULL
  check.lower <- NULL
  
  nComp <- length(object@components)
  if(!is.null(object@calibMoa(normalCopula(0),0))) {
    nonIndep <- sapply(object@components[-nComp], function(x) class(x) != "indepCopula")
    for (i in (1:(nComp-1))[nonIndep]) {
      upParam <- object@calibMoa(object@components[[i]], object@distances[i+1])
      if(any(is.na(upParam))) {
        check.upper <- c(check.upper, TRUE)
      } else {
        if (class(object@components[[i]]) == "frankCopula" && upParam == 0) {
          check.upper <- c(check.upper, TRUE)
        } else {
          check.upper <- c(check.upper, FALSE)
        }
      }
        
      check.lower <- c(check.lower, is.na(object@calibMoa(object@components[[i]], c(0,object@distances)[i])))
    }
    if(sum(check.upper>0)) return(paste("Reconsider the upper boundary conditions of the following copula(s): \n",
                                        paste(sapply(object@components[check.upper], function(x) describeCop(x, "very short")), 
                                              "at", object@distances[check.upper],collapse="\n")))
    if(sum(check.lower>0)) return(paste("Reconsider the lower boundary conditions of the following copula(s): \n",
                                        paste(sapply(object@components[check.lower], function(x) describeCop(x, "very short")), 
                                              "at", object@distances[check.lower],collapse="\n")))
  }
  return(TRUE)
}

setClass("spCopula", representation = representation("copula", 
                                                     components="list",
                                                     distances="numeric", 
                                                     calibMoa="function", 
                                                     unit="character"),
         validity = validSpCopula, contains = list("copula"))

############################
## Spatio-Temporal Copula ##
############################

validStCopula <- function(object) {
  if(length(object@tlags) != length(object@spCopList)) return("The length of the temporal distance vector must equal the number of spatial copulas.")
  return(TRUE) # validity of any spCopula in spCopList is tested by the constructor, I believe
}

setClass("stCopula", representation = representation("copula", 
                                                     spCopList="list", 
                                                     tlags="numeric",
                                                     tres="character"),
         validity = validStCopula, contains = list("copula"))

#########################
## Spatial Vine Copula ##
#########################

validMixedSpVineCopula <- function(object) {
  return(all(sapply(object@spCop,validSpCopula) & validObject(object@topCop)))
}

setClass("mixedSpVineCopula", representation("copula", spCop="list", topCop="copula"),
         validity = validMixedSpVineCopula, contains=list("copula"))

validPureSpVineCopula <- function(object) {
  return(all(sapply(object@spCop,validSpCopula)))
}

setClass("pureSpVineCopula", representation("copula", spCop="list"),
         validity = validPureSpVineCopula, contains=list("copula"))

setClassUnion("spVineCopula",c("mixedSpVineCopula","pureSpVineCopula"))

#################################
## Spatio-temporal Vine Copula ##
#################################

validStVineCopula <- function(object) {
  return(validStCopula(object@stCop) & validObject(object@topCop))
}

setClass("stVineCopula", representation("copula", stCop="stCopula", topCop="copula"),
         validity = validStVineCopula, contains=list("copula"))

########################################
## spatial classes providing the data ##
########################################

## neighbourhood:

validNeighbourhood <- function(object) {
  if(length(var)>1)
    return("Only a single variable name is supported.")
  # check for number of rows
  if (nrow(object@data) != nrow(object@distances)) 
    return("Data and distances have unequal number of rows.")
  if (nrow(object@data) != nrow(object@index)) 
    return("Data and index have unequal number of rows.")
  # check for columns
  if (ncol(object@data) != ncol(object@distances) + 1 + length(object@coVar))
    return("Data and distances have non matching number of columns.")
  if (ncol(object@data) != ncol(object@index) + length(object@coVar)) 
    return("Data and index have non matching number of columns.")
  else 
    return(TRUE)
}

setClass("neighbourhood",
         representation = representation(data = "data.frame", 
                                         distances="matrix", 
                                         index="matrix",
                                         var="character",
                                         coVar="character",
                                         prediction="logical"),         
         validity = validNeighbourhood)

## ST neighbourhood
validStNeighbourhood <- function(object) {
  dimDists <- dim(object@distances)
  dimInd <- dim(object@index)
  
  stopifnot(length(dimDists)==3)
  stopifnot(length(dimInd)==3)
  stopifnot(dimDists[3] == dimInd[3])
  
  if (nrow(object@data) != dimDists[1]) 
    return("Data and distances have unequal number of rows.")
    if (nrow(object@data) != dimInd[1]) 
    return("Data and index have unequal number of rows.")
  if (ncol(object@data) + object@prediction != dimInd[2] + length(object@coVar)) 
    return("Data and index have non matching number of columns.")
  if (dimDists[2]+1 != dimInd[2]) 
    return("Data and index have non matching number of columns.")
  else 
    return(TRUE)
}

setClass("stNeighbourhood",
         representation = representation(data = "data.frame", 
                                         distances="array", 
                                         index="array",
                                         var="character", 
                                         coVar="character",
                                         prediction="logical"),
         validity = validStNeighbourhood)

