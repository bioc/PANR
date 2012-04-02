\name{pvclustModule}
\alias{pvclustModule}
\alias{pvclustModule,PAN,numeric_Or_integer_Or_missing,character_Or_missing,character_Or_missing,logical_Or_missing,logical_Or_missing-method}
\title{
Search enriched functional gene modules by pvclust
}
\description{
The function employs the R package \code{pvclust} to search significant
functional gene modules using hierarchical clustering with bootstrap
resampling. 
}
\usage{
pvclustModule(object, nboot=1000, metric="cosine", hclustMethod="average",
filter=TRUE, verbose=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{nboot}{
a numeric value giving the number of bootstraps for \code{pvclust}.
}
	\item{metric}{
a character value specifying which distance metric to use for the hierarchical
clustering: 'correlation', 'cosine', 'abscor' or those allowed by the argument 
'method' in \code{\link[stats:dist]{dist}}.
}
	\item{hclustMethod}{
the agglomerative method used in hierarchical clustering: 'average', 'ward',
'single', 'complete', 'mcquitty', 'median' or 'centroid' (see the argument
\code{method} in \code{\link[stats:hclust]{hclust}} for more details).
}
	\item{filter}{
a logical value specifying whether or not to filter out screening data of
genes without significant associations with all the other genes.
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
This function performs hierarchical clustering with bootstrap resampling
to quantify significance of gene clusters (modules) based on the package
\code{pvclust}.  
}
\value{
This function will return an object of class \code{PAN} with inferred
gene modules (\code{modules$clusters}) and corresponding p-values (\code{modules$pval})
updated in slot 'modules'.
}
\references{
Xin Wang, Mauro Castro, Klaas W. Mulder and Florian Markowetz, Posterior 
association networks and enriched functional gene modules inferred from 
rich phenotypic perturbation screens, in preparation.

R. Suzuki and H. Shimodaira. Pvclust: an r package for assessing the
uncertainty in hierarchical clustering. Bioinformatics, 22(12):1540, 2006.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
\dontrun{
data(bm, package="PANR")
pan<-new("PAN", bm1=bm1)
pan<-infer(pan, para=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(5)),
filter=FALSE, verbose=TRUE)
data(Bakal2007Cluster, package="PANR")
pan<-buildPAN(pan, engine="igraph", para=list(nodeColor=nodeColor, 
hideNeg=TRUE), verbose=TRUE)
##need pvclust to search modules
library(pvclust)
pan<-pvclustModule(pan, nboot=1000, metric="cosine",
hclustMethod="average", filter=TRUE, verbose=TRUE, r=c(5:12/7))
}
}












