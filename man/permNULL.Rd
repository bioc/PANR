\name{permNULL}
\alias{permNULL}
\alias{permNULL,BetaMixture,character_Or_missing-method}
\title{
Do permutations for input rich phenotyping screens.
}
\description{
This is an internal S4 method for class \code{PAN} to permute the original 
inputted phenotyping screens to further fit the dissociation component of 
the beta-mixture model. 
}
\usage{
permNULL(object, permMethod="keepRep", ...)
}
\arguments{
	\item{object}{
	an object of \code{BetaMixture}.
}
	\item{permMethod}{
	a character value of the method to permute, either 'keepRep' (keep
	the order of replicates) or 'all' (regardless of replicates).
}
}
\value{
a matrix of permuted data.
}
\seealso{
\code{\link[PANR:fitNULL]{fitNULL}}, \code{\link[PANR:fitBM]{fitBM}}, \code{\link[PANR:BetaMixture]{BetaMixture}}
}
\references{
Xin Wang, Roland F. Schwarz, Mauro Castro, Klaas W. Mulder 
and Florian Markowetz, Posterior association networks and enriched 
functional gene modules inferred from rich phenotypic perturbation 
screens, in preparation.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
data(bm)
bm.perm<-permNULL(bm1)
}
\keyword{internal}












