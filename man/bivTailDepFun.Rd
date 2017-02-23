\name{bivJointDepFun}
\alias{bivJointDepFun}
\alias{lowerBivJointDepFun}
\alias{upperBivJointDepFun}

\alias{empBivJointDepFun}
\alias{lowerEmpBivJointDepFun}
\alias{upperEmpBivJointDepFun}

\title{
Bivariate joint dependence functions
}
\description{
Functions returning a (empirical) joint shortfall/exceedance dependence function for a bivariate copula (sample). The functions can be lower (shortfall), upper (exceedance), or both, where values below 0.5 are calculated from the lower joint dependence function and values larger 0.5 from the upper joint dependence function. The definition follows the one by Nelsen for the tail dependence index (2006, Theorem 5.4.2.):

lower joint dependence (shortfall):
\eqn{\frac{P(U \leq u | V \leq u)}{P(V \leq v)}=\frac{C(u,u)}{u}}{P(U <= u | V <= u)/P(V <= u) = C(u,u)/u}

upper joint dependence (exceedance):
\eqn{\frac{P(U \geq u | V \geq u)}{P(V \geq v)}=\frac{1-2u+C(u,u)}{1-u}}{P(U >= u | V >= u)/P(V >= u) = (1 - 2u + C(u,u))/(1-u)}

Taking the limits of the above for \eqn{u \rightarrow 0}{u -> 0} and \eqn{u \rightarrow 1}{u -> 1} respectively yields the lower and upper tail dependence index. For tail dependence of multivariate copulas, see the papers by Joe et al. (2010) and Nikoloulopoulos et al. (2012).
}
\usage{
bivJointDepFun(copula)
lowerBivJointDepFun(copula)
upperBivJointDepFun(copula)

empBivJointDepFun(u)
lowerEmpBivJointDepFun(u)
upperEmpBivJointDepFun(u)
}

\arguments{
  \item{copula}{
an object of class \code{\linkS4class{copula}}
}
  \item{u}{
a bivariate sample on (0,1)
}
}

\value{
A function taking arguments from the unit interval (0,1) and returning the corresponding index.
}
\references{
Definition: 

Nelsen, Roger B. An introduction to copulas. Vol. 139. Springer Science & Business Media, 2013.

Further reading on multivariate tail dependence:

Joe, Harry, Haijun Li, and Aristidis K. Nikoloulopoulos. "Tail dependence functions and vine copulas." Journal of Multivariate Analysis 101.1 (2010): 252-270.

Nikoloulopoulos, Aristidis K., Harry Joe, and Haijun Li. "Vine copulas with asymmetric tail dependence and applications to financial return data." Computational Statistics & Data Analysis 56.11 (2012): 3659-3673.

Plots inspired by:

\url{http://freakonometrics.hypotheses.org/2435}
}
\author{
Benedikt Graeler
}

\examples{
library("VineCopula")
data("simulatedTriples")
X <- rankTransform(triples[,c(1,3)])
  
tdfEmp <- empBivJointDepFun(X)
plot(tdfEmp,ylim=c(0,1),
     ylab="tail dependence index")  
abline(v=0.5, col="grey")

smplTau <- cor(X,method="kendall")[1,2]

# Gauss
tdfGauss <- bivJointDepFun(normalCopula(sin(smplTau*pi/2)))
curve(tdfGauss,add=TRUE,col="blue")

# survival Gumbel
tdfGumbel <- bivJointDepFun(surGumbelCopula(1/(1-smplTau)))
curve(tdfGumbel,add=TRUE,col="darkgreen")

# survival BB6 copula
tdfBB6 <- bivJointDepFun(surBB6Copula(c(4.65,2.28)))
curve(tdfBB6,add=TRUE,col="red")

legend("bottomleft",c("empircal","Gauss","surv. Gumbel","surv. BB6"),
       col=c("black","blue","darkgreen","red"),lty=1)
} 