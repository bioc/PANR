\name{BetaMixture-class}
\alias{BetaMixture}
\alias{BetaMixture-class}
\docType{class}
\title{
An S4 class for beta mixture modelling of functional gene associations
}
\description{
This S4 class includes methods to do beta-mixture modelling of functional
gene associations given rich phenotyping screens. 
}
\section{Objects from the Class}{ 
Objects of class \code{BetaMixture} can be created from \code{new("BetaMixture",
metric, order, association, model, pheno, partition)} (see the example
below for details).
}
\section{Slots}{
	\describe{	
		\item{\code{pheno}:}{
		a numeric matrix of rich phenotypes with rows and columns
		specifying genes and samples, respectively.
		}
		\item{\code{metric}:}{
		a character value specifying the metric to compute similarity
		scores. Currently, 'cosine' and 'correlation' are supported (see
		\code{\link[PANR:assoScore]{assoScore}} for more details).
		}
		\item{\code{order}:}{
		a numeric value specifying the order of the similarity score to
		be computed. Only 1 and 2 is supported for the current version.
		The first order (when \code{order=1}) similarities are used for
		quatification of the strength of functional associations between
		genes, whilst the second order (when \code{code=2}) ones are employed
		to compute the strength of modularity between genes.
		}
		\item{\code{association}:}{
		a numeric vector providing all association scores between genes.
		This can be useful when \code{pheno} is not available or the user
		has a different way to compute functional associations.
		}
		\item{\code{model}:}{
		a character value specifying whether the original (if \code{global})
		or extended (if \code{stratified}) model is used. 
		}
		\item{\code{partition}:}{
		a numeric of gene partition labels (e.g. c(rep(1, 100), rep(2, 20), 
		rep(3, 80)) is a valid vector of partition labels for a vector of 
		associations falling into three categories of interaction types 1, 
		2 and 3).
		}
		\item{\code{result}:}{
		a list storing results from S4 methods of this class.
		}
		\item{\code{summary}:}{
		a list of summary information for available results.
		}
	}
}
\section{Methods}{
	An overview of methods (More detailed introduction can be found in
	help for each specific function.):
	\describe{	
		\item{\code{permNULL}}{
		do permutations for input rich phenotyping screens (\code{'pheno'}).
		}
		\item{\code{fitNULL}}{
		fit the permuted association scores to a beta distribution. 
		}
		\item{\code{fitBM}}{
		fit the functional association scores computed from input screens
		to a three-beta mixture model.
		}
		\item{\code{view}}{
		view the fitting results (a histogram of the original data and
		fitted probability density curves) for NULL and real data.
		}
		\item{\code{summarize}}{
		summarize results including input data and parameters, NULL fitting
		and beta mixture fitting.
		}
	}
}
\author{Xin Wang \email{xw264@cam.ac.uk}}
\references{
Xin Wang, Roland F. Schwarz, Mauro Castro, Klaas W. Mulder 
and Florian Markowetz, Posterior association networks and enriched 
functional gene modules inferred from rich phenotypic perturbation 
screens, in preparation.
}
\seealso{
  \code{\link[PANR:permNULL]{permNULL}}
  \code{\link[PANR:fitNULL]{fitNULL}}
  \code{\link[PANR:fitBM]{fitBM}}
  \code{\link[PANR:view]{view}}
  \code{\link[PANR:summarize]{summarize}}
}
\examples{
\dontrun{
	data(Bakal2007)
	bm1<-new("BetaMixture", pheno=Bakal2007, metric="cosine",
		model="global", order=1)
	bm1<-fitNULL(bm1, nPerm=10, thetaNULL=c(alphaNULL=4, betaNULL=4),
		sumMethod="median", permMethod="all", verbose=TRUE)
	bm1<-fitBM(bm1, para=list(zInit=NULL, thetaInit=c(alphaNeg=2, betaNeg=4, 
		alphaNULL=bm1@result$fitNULL$thetaNULL[["alphaNULL"]], 
		betaNULL=bm1@result$fitNULL$thetaNULL[["betaNULL"]], 
		alphaPos=4, betaPos=2), gamma=NULL), 
		ctrl=list(fitNULL=FALSE, tol=1e-1), verbose=TRUE, gradtol=1e-3)
	view(bm1, "fitNULL")
	view(bm1, "fitBM")
	bm1
}
}
\keyword{classes}
