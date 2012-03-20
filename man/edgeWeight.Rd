\name{edgeWeight}
\alias{edgeWeight}
\alias{edgeWeight, PAN-method}
\title{
Compute edge weights for posterior association networks
}
\description{
This is an internal function to compute edge weights before 
inferring a posterior association network. 
}
\usage{
edgeWeight(object, which="bm1", type="SNR", log=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{which}{
a character value specifying which BetaMixture modelling result to use:
first-order (if 'bm1') or second-order (if 'bm2').
}
	\item{type}{
a character value giving the type of edge weight to compute: signal-
to-noise ratio (if 'SNR'), posterior probability odd (if 'PPR') or
posterior probability (if 'PP').
}
	\item{log}{
a logical value specifying whether or not to compute logrithms for edge
weights.
}
}
\details{
This function will be called by \code{\link[PAN:infer]{infer}} to compute edge weights for posterior
association networks. When inferring a signed PAN, signal-to-noise ratios
are suggested to use; while inferring a PAN of only positive associations,
posterior probability odds or posterior probabilities are preferred.
}
\value{
This function will return a numeric adjacency matrix of edge weights.
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
\seealso{
  \code{\link[PAN:infer]{infer}}
}

\keyword{internal}











