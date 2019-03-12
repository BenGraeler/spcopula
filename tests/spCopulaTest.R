library("spcopula")
data("spCopDemo")
calcKTauPol <- fitCorFun(bins, degree=3)

spCop <- spCopula(components=list(normalCopula(), tCopula(),
                                  frankCopula(), normalCopula(), claytonCopula(),
                                  claytonCopula(), claytonCopula(), claytonCopula(),
                                  claytonCopula(), indepCopula()),
                  distances=c(0, bins$meanDists[1:9]),
                  spDepFun=calcKTauPol, unit="m")


invdduCopula(0.5, spCop, 0.3, h=300)
invddvCopula(0.5, spCop, 0.3, h=300)

y <- dduCopula(c(0.3,0.7), spCop, h=300)
invdduCopula(0.3, spCop, y, h=300)

y <- ddvCopula(c(0.3,0.7), spCop, h=300)
invddvCopula(0.7, spCop, y, h=300)


