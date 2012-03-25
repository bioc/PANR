\name{fitBM}
\alias{fitBM}
\alias{fitBM,BetaMixture,list_Or_missing,list_Or_missing,logical_Or_missing-method}
\title{
Fit a three-beta mixture model to densities of functional gene associations
}
\description{
The function fits a three-beta mixture model to densities of functional
gene associations computed from rich phenotyping screens.
}
\usage{
fitBM(object, para=list(zInit=NULL, thetaInit=c(alphaNeg=2, betaNeg=4,
alphaNULL=4, betaNULL=4, alphaPos=4, betaPos=2), gamma=NULL), 
ctrl=list(fitNULL=FALSE, tol=1e-3, maxIter=NULL), verbose=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{BetaMixture}.
}
	\item{para}{
a list of initial values for parameter estimation in fitting a three-beta
mixture model (see 'details'). 
}
	\item{ctrl}{
a list of control parameters for the mixture model fitting (see 'details').
}
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
other arguments of the function \code{\link[stats:nlm]{nlm}}.
}
}
\details{
This function fits a beta-mixture model to functional gene associations
using the Expectation-Maximization algorithm. The function allows various
parameter settings to perform fitting by the original (if model='global')
or stratified (if model='stratified') beta-mixture model (the model should
be specified when creating a new object of \code{\link[PANR:BetaMixture]{BetaMixture}}).

The initial values of the beta distributions can be set by \code{thetaInit},
is a numeric vector including the two shape parameters for the '-'
(negative), 'x' (NULL) and '+' (positive) distributions. Please note that
if \code{ctrl$NULL} is \code{TRUE}, meaning that the NULL distribution
has already been fitted, then \code{para$alphaNULL} and \code{para$betaNULL}
are supposed to be filled in the estimated NULL parameters by the function
\code{\link[PANR:fitNULL]{fitNULL}}).

\code{zInit} is a matrix of posterior probabilities for gene associations
following the three mixture components. 

The hyper-parameters for the dirichlet priors for the mixture components
can also be set by \code{para$gamma}, which is a numeric matrix with rows
and columns correponding to association partitions and the three beta
mixture components.

The other arguments to control the fitting algorithm are \code{tol} and \code{maxIter},
which are covergence tolerence and the maximal iterations.

Since the estimation of shape parameters of beta distributions are realized
by the function \code{\link[stats:nlm]{nlm}} numerically, additional arguments for \code{\link[stats:nlm]{nlm}}
are allowed by \code{...}.
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
  \code{\link[PANR:fitNULL]{fitNULL}}
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
}
}
