library("spcopula")
spCopT0 <- spCopula(components=list(claytonCopula(8), claytonCopula(4), 
                                    claytonCopula(2), claytonCopula(1),
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300,400,500),
                    unit="km")
spCopT1 <- spCopula(components=list(claytonCopula(4), claytonCopula(2), 
                                    claytonCopula(1), claytonCopula(0.5),
                                    indepCopula()),
                    distances=c(0,100,200,300,400),
                    unit="km")
spCopT2 <- spCopula(components=list(claytonCopula(2), claytonCopula(1), 
                                    claytonCopula(0.5), indepCopula()),
                    distances=c(0,100,200,300),
                    unit="km")

stCop <- stCopula(components=list(spCopT0, spCopT1, spCopT2),
                  tlags=-(0:2))

invdduCopula(0.5, stCop, 0.3, h=matrix(c(150,-1),ncol=2))
invddvCopula(0.5, stCop, 0.3, h=matrix(c(150,-1),ncol=2))

y <- dduCopula(c(0.3,0.7), stCop, h=matrix(c(150,-1),ncol=2))
invdduCopula(0.3, stCop, y, h=matrix(c(150,-1),ncol=2))

y <- ddvCopula(c(0.3,0.7), stCop, h=matrix(c(150,-1),ncol=2))
invddvCopula(0.7, stCop, y, h=matrix(c(150,-1),ncol=2))

