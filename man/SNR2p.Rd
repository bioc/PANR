\name{SNR2p}
\alias{SNR2p}
\alias{SNR2p,BetaMixture,numeric_Or_integer-method}
\title{
Translate p-values to Signal-to-Noise Ratios
}
\description{
The function translate Signal-to-Noise Ratios to p-values based on the fitted 
mixture model.
}
\usage{
SNR2p(object, SNR, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{BetaMixture}.
}
	\item{SNR}{
a numeric or integer specifying the SNR to translate
}
	\item{...}{
not in use, only for further extension.
}
}
\value{
The function will return a data frame including SNRs, lower and upper 
quantiles and corresponding p-values.
}
\references{
Xin Wang, Mauro Castro, Klaas W. Mulder and Florian Markowetz, Posterior 
association networks and enriched functional gene modules inferred from 
rich phenotypic perturbation screens, in preparation.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\seealso{
  \code{\link[PANR:SNR2p]{SNR2p}}
}
\examples{
\dontrun{
data(Bakal2007)
bm1<-new("BetaMixture", pheno=Bakal2007, model="global", order=1)
bm1<-fitNULL(bm1, nPerm=10, thetaNULL=c(alphaNULL=4, betaNULL=4),
	sumMethod="median", permMethod="all", verbose=TRUE)
bm1<-fitBM(bm1, para=list(zInit=NULL, thetaInit=c(alphaNeg=2, betaNeg=4, 
	alphaNULL=bm1@result$fitNULL$thetaNULL[["alphaNULL"]], 
	betaNULL=bm1@result$fitNULL$thetaNULL[["betaNULL"]], 
	alphaPos=4, betaPos=2), gamma=NULL), 
	ctrl=list(fitNULL=FALSE, tol=1e-1), verbose=TRUE, gradtol=1e-3)
SNR2p(bm1, SNR=10)
}
}
