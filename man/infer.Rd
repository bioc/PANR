\name{infer}
\alias{infer}
\alias{infer,PAN,list,logical_Or_missing,logical_Or_missing-method}
\title{
Infer a posterior association network
}
\description{
The function infers a posterior association network from beta-mixture
modelling of functional associations computed from rich phenotyping screens. 
}
\usage{
infer(object, para=list(type='SNR', log=TRUE, sign=TRUE, cutoff=0),
filter=FALSE, verbose=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{para}{
a list of parameters to perform inference (see details).
}
	\item{filter}{
a logical value specifying whether or not to filter out genes without any
significant association with all the other genes.
}
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
not in use, only for further extension.
}
}
\details{
This function employs different edge weights to infer a posterior association
network (see \code{\link[PAN:edgeWeight]{edgeWeight}} for more details).
Multiple parameters are provided for the user to specify the network:

'type' -  a character value giving the type of edge weights: signal-to-
noise ratio ('SNR'), posterior probability ratio ('PPR') or posterior
probability ('PP').

'log' - a logical value specifying whether or not to compute logrithms
for edge weights.

'sign' - a logical value specifying whether a signed graph should be
inferred. It is only valid when \code{type='SNR'}.

'cutoff' -  a numeric value giving the threshold to tell the significance
of an edge.
}
\value{
This function will return an object of class \code{PAN} with inferred PAN
updated in slot 'graph'.
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
  \code{\link[PAN:edgeWeight]{edgeWeight}}
}

\examples{
data(bm, package="PANR")
pan<-new("PAN", bm1=bm1)
pan<-infer(pan, para=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(5)),
filter=FALSE, verbose=TRUE)
}












