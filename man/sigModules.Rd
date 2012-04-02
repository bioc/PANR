\name{sigModules}
\alias{sigModules}
\alias{sigModules,PAN,numeric_Or_integer_Or_missing,numeric_Or_integer_Or_missing,numeric_Or_integer_Or_missing,character_Or_missing,logical_Or_missing-method}
\title{
Retrieve ids for significant gene modules searched by pvclust
}
\description{
The function retrieve ids for significant gene modules that satisfy the 
given p-value cutoff and module size range. 
}
\usage{
sigModules(object, pValCutoff=0.01, minSize=3, maxSize=100, sortby="size", 
decreasing=FALSE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{pValCutoff}{
a numeric value specifying the p-value cutoff to tell the significance of a 
gene module.
}
	\item{minSize}{
a numeric or integer value giving the minimal size of gene modules.
}
	\item{maxSize}{
a numeric or integer value giving the maximal size of gene modules.
}
	\item{sortby}{
a character value specitying how to sort the list of gene modules: by `size' 
(module size) or `pval' (pvclust p-value).
}
	\item{decreasing}{
a logical value specifying whether or not the sorting will be conducted 
decreasingly.
}
	\item{...}{
not in use, but only for further extension.
}
}
\details{
This function facilitates the user to retrieve significant gene modules found 
by pvclust and obtain their ids, which can be subsequently used for visualization 
(see \code{\link[PANR:viewPAN]{viewPAN}} for details).  
}
\value{
a numeric vector of ids for significant gene modules
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
inds<-sigModules(pan, pValCutoff=0.01, minSize=5, maxSize=100, sortby="size", 
decreasing=FALSE)
pan@modules$clusters[inds]
}
}












