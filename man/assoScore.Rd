\name{assoScore}
\alias{assoScore}
\title{
Association scores for gene pairs
}                                                                        
\description{
This function compute similarity scores to quantify associations between
pairs of genes, given measured rich phenotypes.
}
\usage{
assoScore(pheno, metric="cosine", upperTri=TRUE, transform=TRUE,
verbose=TRUE, ...) 
}
\arguments{
	\item{pheno}{
a numeric matrix of rich phenotypes with rows and columns specifying
samples and genes, respectively.
}
	\item{metric}{
a character value specifying the metric to compute similarity scores.
Currently, 'cosine' and 'correlation' are supported (see details for
more).
}
	\item{upperTri}{
a logical value specifying whether (if \code{TRUE}) to take the upper
triangular of the similarity matrix or not (if \code{FALSE}).
}
	\item{transform}{
a logical value specifying whether to transform (if \code{TRUE}) the range of
association scores from \code{[-1, 1]} to \code{[0, 1]} or not (if \code{FALSE}).
}
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
other arguments for function \code{cor}
}
}
\details{
This function aims at quantifying the associations between genes of
interest given certain phenotyping measurements (e.g. gene expressions
by microarray, cell viabilities, morphological phenotypes, etc.). For
the current version of the package, the user can either choose
'cosine' or 'correlation'. When the latter is chosen, additional
arguments (e.g. 'method') for the function \code{cor} are allowed. 
}
\value{
This function will return either a vector (if \code{upperTri=TRUE}) or a matrix
(if \code{upperTri=FALSE}) of association scores for given phenotypes. 
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
\code{\link[stats:cor]{cor}}, \code{\link[PAN:cosineSim]{cosineSim}}
}
\examples{
toydata<-matrix(rnorm(n=2000, mean=0, sd=4), nrow=100, ncol=20)
toyasso<-assoScore(t(toydata), "cosine", upperTri=FALSE, transform=FALSE)
##transform to [0, 1]
toyasso01<-assoScore(t(toydata), "cosine", upperTri=FALSE, transform=TRUE)
##transform to [0, 1] and return only the upper triangular
toyasso01upper<-assoScore(t(toydata), "cosine", upperTri=TRUE, transform=
TRUE)
##use spearman correlation
toyassoSp<-assoScore(t(toydata), "correlation", upperTri=FALSE, transform=
FALSE, method="spearman")
}










