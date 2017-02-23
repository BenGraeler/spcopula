## stVineCopula estimation following the vignette, but using a randomly 
## selected smaller subset of the original data to reduce calculation demands.
## Thus, results are likely to differ (a little) from the original study.

library("VineCopula")
library("spcopula")
data("EU_RB_2005")

## spatio-temporal copula
# binning, using only 90 out of all temporal instances
stBins <- calcBins(EU_RB_2005, "rtPM10", nbins=40, t.lags=-(0:2),
                   instances=10)

calcKTau <- fitCorFun(stBins,c(3,3,3))

families <- c(normalCopula(0), tCopula(0), claytonCopula(0), 
              frankCopula(1), gumbelCopula(1), joeBiCopula(), 
              cqsCopula(), asCopula())

loglikTau <- list()
for(j in 1:3) { # j <-1
  tmpBins <- list(meanDists=stBins$meanDists,
                  lagData=lapply(stBins$lagData, function(x) x[,c(2*j-1,2*j)]))
  loglikTau[[j]] <- loglikByCopulasLags(tmpBins, families, calcKTau[[j]])
}

bestFitTau <- list()
bestFitTau[[1]] <- apply(apply(loglikTau[[1]]$loglik[,-8], 1, rank),
                         2, which.max)
bestFitTau[[2]] <- apply(apply(loglikTau[[2]]$loglik, 1, rank),
                         2, which.max)
bestFitTau[[3]] <- apply(apply(loglikTau[[3]]$loglik, 1, rank),
                         2, which.max)

# gather the fitted copulas and representative distances
listCops <- NULL

for(t.level in 1:3) {
  listCops[[t.level]] <- sapply(1:35, 
                                function(i) {
                                  tmpTau <- calcKTau[[t.level]](stBins$meanDists[i])
                                  if (tmpTau < 0.05) {
                                    return(indepCopula(dim=2))
                                  } else {
                                    return(loglikTau[[t.level]]$copulas[[bestFitTau[[t.level]][i]]][[i]])
                                  }})
}

listDists <- NULL
listDists[[1]] <- stBins$meanDists[1:35]
listDists[[2]] <- stBins$meanDists[1:35]
listDists[[3]] <- stBins$meanDists[1:35]

# build the spatio-temporal copula
stConvCop <- stCopula(components=listCops, distances=listDists,
                      t.lags=c(0,-1,-2))

# get the neighbours
stNeigh <- getStNeighbours(EU_RB_2005, var="rtPM10", spSize=4, 
                           t.lags=-(0:2), timeSteps=10, min.dist=10)


stVineFit <- fitCopula(stVineCopula(stConvCop,vineCopula(9L)), stNeigh,
                       method="indeptest")
stVine <- stVineFit@copula

# retrieve the log-likelihood
stVineFit@loglik