\name{cosineSim}
\alias{cosineSim}
\alias{cosineDist}
\title{
Compute cosine similarities or distances between pairs of genes
}                                                                        
\description{
This function compute cosine similarities or distances between pairs of
genes, given measured rich phenotypes.
}
\usage{
cosineSim(x)
cosineDist(x)
}
\arguments{
	\item{x}{
a numeric matrix of rich phenotypes with rows and columns specifying
samples and genes, respectively.
}
}
\value{
This function will return either a numeric matrix of cosine similarities
or an object of \code{dist}. 
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
\code{\link[PAN:assoScore]{assoScore}}
}
\examples{
toydata<-matrix(rnorm(n=2000, mean=0, sd=4), nrow=100, ncol=20)
toycossim<-cosineSim(t(toydata))
toycosdist<-cosineDist(t(toydata))
}










