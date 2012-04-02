\name{p2SNR}
\alias{p2SNR}
\alias{p2SNR,BetaMixture,numeric_Or_integer}
\title{
Translate p-values to Signal-to-Noise Ratios
}
\description{
The function translate p-values to Signal-to-Noise Ratios based on the fitted 
mixture model.
}
\usage{
p2SNR(object, pval, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{BetaMixture}.
}
	\item{pval}{
a numeric or integer specifying the p-value to translate
}
	\item{...}{
not in use, only for further extension.
}
}
\value{
The function will return a data frame including p-values, lower and upper 
quantiles and corresponding Signal-to-Noise Ratios.
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
  \code{\link[PANR:p2SNR]{p2SNR}}
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
p2SNR(bm1, pval=0.01)
}
}
