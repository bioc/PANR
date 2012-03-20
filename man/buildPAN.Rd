\name{buildPAN}
\alias{buildPAN}
\alias{buildPAN,PAN,character_Or_missing,list_Or_missing,logical_Or_missing-method}
\title{
Build an igraph or RedeR graph for PAN
}
\description{
The function builds a graph for the inferred PAN so that it can be visualize 
in \code{\link[igraph:igraph]{igraph}} or \code{\link[RedeR:RedeR]{RedeR}}
}
\usage{
buildPAN(object, engine="igraph", para=list(nodeColor=NULL, nodeSize=NULL, 
edgeWidth=NULL, edgeColor=NULL, nodeSumCols=1, nodeSumMethod="none", 
hideNeg=TRUE), verbose=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{engine}{
a character value specifying which graphics engine to use: `igraph' or `RedeR'.
}
	\item{para}{
a list of parameters specifying graph attributes (see details)
}
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
all the other arguments accepted by the function \code{\link[pvclust:pvclust]{pvclust}}.
}
}
\details{
Here are the introductions for the detailed graph attributes that the user 
can specify:

'nodeColor' and 'nodeSize' - a vector of node colors or sizes. Please note
that the order of color or size must be concordant with the gene ids in
the argument \code{pheno} of the object of \code{\link[PAN:BetaMixture]{BetaMixture}}.

'nodeSumCols' and 'nodeSumMethod' - these two arguments are used to scale
the colors of nodes by phenotypes. The former argument is a numeric vector
specifying the columns in slot \code{pheno} of the object of class \code{\link[PAN:BetaMixture]{BetaMixture}};
while the latter one is a character value giving the method to summarize
these columns of phenotypes: either 'mean' or 'median'.

'edgeColor' and 'edgeWidth' - a vector of edge colors or width.

'pValCutoff' - the argument is only used when \code{what='module'} and
for module searching based on \code{\link[pvclust:pvclust]{pvclust}}. Only
significant modules will be displayed.

'minSize' and 'maxSize' - two arguments controlling the size of modules
which will be used to filtered out modules that are too small or too large.

'hideNeg' - a logical value specifying whether or not to hide edges with
negative associations
}
\value{
This function will return an object of class \code{PAN} with inferred
gene modules (\code{modules$clusters}) and corresponding p-values (\code{modules$pval})
updated in slot 'modules'.
}
\references{
Xin Wang, Roland F. Schwarz, Mauro Castro, Klaas W. Mulder 
and Florian Markowetz, Posterior association networks and enriched 
functional gene modules inferred from rich phenotypic perturbation 
screens, in preparation.

R. Suzuki and H. Shimodaira. Pvclust: an r package for assessing the
uncertainty in hierarchical clustering. Bioinformatics, 22(12):1540, 2006.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
data(bm, package="PANR")
pan<-new("PAN", bm1=bm1)
pan<-infer(pan, para=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(5)),
filter=FALSE, verbose=TRUE)
data(Bakal2007Cluster, package="PANR")
pan<-buildPAN(pan, engine="igraph", para=list(nodeColor=nodeColor, 
hideNeg=TRUE), verbose=TRUE)
}












