\name{summarize}
\alias{summarize}
\alias{summarize,BetaMixture,character_Or_missing-method}
\alias{summarize,PAN,character_Or_missing-method}
\title{
Summarize the object of S4 class 'BetaMixture' or 'PAN'
}
\description{
The function helps print a summary of an object of S4 class \code{BetaMixture}
or \code{PAN}.
}
\usage{
summarize(object, what='ALL', ...) 
}
\arguments{
	\item{object}{
an object of S4 class \code{BetaMixture} or \code{PAN}.
}
	\item{what}{
a character value specifying what to print (see details).
}
	\item{...}{
not in use, only for further extension.
}
}
\details{
This function print a summary of an object of \code{BetaMixture} or \code{PAN}.
The function is also called by S4 method \code{show}, which prints only
a short message about the input parameters and data.

For an object of class \code{BetaMixture}:

If \code{what='input'}, the function prints to screen a summary of input
parameters;
If \code{what='fitNULL'}, the function prints to screen a summary of
fitting results for the NULL distribution.
If \code{what='fitBM'}, the function prints to screen a summary of
fitting results for the beta-mixture model.
If \code{what='ALL'}, all above messages will be printed.

For an object of class \code{PAN}:

If \code{what='input'}, the function prints to screen a summary of input
object(s) of class \code{BetaMixture};
If \code{what='graph'}, the function prints to screen a summary of inferred
posterior association network;
If \code{what='module'}, the function prints to screen a summary of 
functional gene modules;
If \code{what='ALL'}, all above messages will be printed.
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

\examples{
data(bm)
summarize(bm1, what='ALL')
}












