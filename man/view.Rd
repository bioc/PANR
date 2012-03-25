\name{view}
\alias{view}
\alias{view,BetaMixture,character_Or_missing-method}
\title{
View the results of beta-mixture model fitting
}
\description{
This function facillitate the user to view and check the fitting results
of NULL and beta-mixture distributions.
}
\usage{
view(object, what="fitBM", ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{BetaMixture}.
}
	\item{what}{
a character value specifying to show the fitting results of NULL (if
\code{what='fitNULL'}) or beta-mixture model (if \code{what='fitBM'}).
}
	\item{...}{
not in use, only for further extension.
}
}
\details{
The function help the user to view and check the fitting of NULL and
three beta-mixture model to permuted and real phenotyping screens,
respectively. For either fitting, a histogram of association scores
and the fitted beta distribution(s) will be plotted. For beta-mixture
fitting, the integrated probability density function for the mixed
distribution will also be plotted. 
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}
\references{
Xin Wang, Roland F. Schwarz, Mauro Castro, Klaas W. Mulder 
and Florian Markowetz, Posterior association networks and enriched 
functional gene modules inferred from rich phenotypic perturbation 
screens, in preparation.
}
\seealso{
  \code{\link[PANR:fitNULL]{fitNULL}}
  \code{\link[PANR:fitNULL]{fitBM}}
}

\examples{
data(bm)
view(bm1, what="fitNULL")
view(bm1, what="fitBM")
}












