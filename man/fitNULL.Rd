\name{fitNULL}
\alias{fitNULL}
\alias{fitNULL,BetaMixture,numeric_Or_integer_Or_missing,numeric_Or_integer_Or_missing,character_Or_missing,character_Or_missing,logical_Or_missing-method}
\title{
Fit the NULL component of a three-beta mixture model for functional gene
associations
}
\description{
The function performs permutations to the input rich phenotyping screens,
and subsequently fit a beta distribution to the densities. 
}
\usage{
fitNULL(object, nPerm=20, thetaNULL=c(alphaNULL=4, betaNULL=4),
sumMethod="median", permMethod="keepRep", verbose=TRUE, ...) 
}
\arguments{
	\item{object}{
an object of S4 class \code{BetaMixture}.
}
	\item{nPerm}{
a positive numeric or integer value specifying the number of permutations.
}
	\item{thetaNULL}{
a list of numeric values giving the initial values for estimating
the two shape parameters (see \code{\link[stats:dbeta]{dbeta}} for more details).
}
	\item{sumMethod}{
a character value specifying how to summarize estimated parameters from
multiple permutations. The current version only supports 'median' and
'mean'.
}
	\item{permMethod}{
a character value of the method to permute, either 'keepRep' (keep the
order of replicates) or 'all' (regardless of replicates).
}
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
other arguments for function \code{\link[MASS:fitdistr]{fitdistr}}.
}
}
\details{
This function is a step prior to the three-beta mixture model fitting to
functional gene association scores. The fitted parameters are then used
as a fixed parameters in the three-beta mixture model for further fitting
to the real screens. The NULL fitting is performed using the function
\code{\link[MASS:fitdistr]{fitdistr}}, so other arguments for the function fitdistr are also allowed
by the argument \code{...}.
}
\value{
This function will return an updated object of class \code{BetaMixture}.
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
  \code{\link[PAN:fitBM]{fitBM}}
}

\examples{
data(Bakal2007)
bm1<-new("BetaMixture", pheno=Bakal2007, model="global", order=1)
bm1<-fitNULL(bm1, nPerm=10, thetaNULL=c(alphaNULL=4, betaNULL=4),
sumMethod="median", permMethod="all", verbose=TRUE)
}












