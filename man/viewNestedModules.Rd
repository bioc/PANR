\name{viewNestedModules}
\alias{viewNestedModules}
\alias{viewNestedModules,PAN,numeric_Or_integer_Or_missing,numeric_Or_integer_Or_missing,numeric_Or_integer_Or_missing,logical_Or_missing-method}
\title{
View the nested modules in a posterior association network in RedeR
}
\description{
The function displays the nested enriched functional gene modules found by 
pvclust in a powerful graphic visualization software \code{\link[RedeR:RedeR]{RedeR}}.
}
\usage{
viewNestedModules(object, pValCutoff=0.01, minSize=3, maxSize=100, 
verbose=TRUE, ...)
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
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
not in use, but only for further extension.
}
}
\details{
This function presents the searched enriched functional modules in \code{\link[RedeR:RedeR]{RedeR}} 
- a bioconductor package for network visualization.

Please note that the user is expected to run \code{\link[PANR:buildPAN]{buildPAN}} to build a 
graph and search modules using \code{\link[PANR:pvclustModule]{pvclustModule}} prior to visualize using 
this function.

Please also note that if `RedeR' is selected as the graphics engine, it is 
suggested to manually organise the sizes and positions of containers (for 
nesting gene modules) run a dynamic layout to obtain the best structure for 
the network. 
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
  \code{\link[RedeR:addGraph]{addGraph}}, \code{\link[RedeR:nestNodes]{nestNodes}}, \code{\link[PANR:viewPAN]{viewPAN}}, \code{\link[PANR:buildPAN]{buildPAN}}
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
pan<-pvclustModule(pan, nboot=10000, metric="cosine2",
hclustMethod="average", filter=TRUE, verbose=TRUE, r=c(5:12/7))
viewNestedModules(pan, pValCutoff=0.05, minSize=5, maxSize=100)
}
}













