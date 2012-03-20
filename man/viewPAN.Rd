\name{viewPAN}
\alias{viewPAN}
\alias{viewPAN,PAN,character_Or_missing,numeric_Or_integer_Or_missing,character_Or_missing,logical_Or_missing-method}
\title{
Show posterior association networks or modules in `igraph' or `RedeR'
}
\description{
The function display inferred posterior association networks or enriched
functional gene modules in igraph or a powerful graphical visualization
software \code{\link[RedeR:RedeR]{RedeR}}.
}
\usage{
viewPAN(object, what="graph", moduleID=1,
layout="layout.fruchterman.reingold", verbose=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{what}{
a character value specifying which to show: 'graph' or 'pvclustModule'.
}
	\item{moduleID}{
a numeric or integer vector of modules to view (see details in \code{\link[PAN:sigModules]{sigModules}}).
This argument will be applied only when what=`pvclustModule'.
}
	\item{layout}{
a character value specifying the layout method (see details in \code{\link[igraph:layout]{layout}}).
This argument will be applied only when engine=`igraph' when build a graph
for PAN using \code{\link[PAN:buildPAN]{buildPAN}}.
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
This function presents the inferred posterior association network or enriched 
functional modules in igraph or \code{\link[RedeR:RedeR]{RedeR}} depending on the graphics
engine used when building the graph for PAN (details in \code{\link[PAN:buildPAN]{buildPAN}}).

Please note that when viewing a dense PAN in `igraph', it could be very
messy. Multiple modules can be viewed at the same time in `RedeR' powered
by its feature of `containers', which are used to group gene modules. When
view multiple modules in `igraph', 

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
  \code{\link[PAN:buildPAN]{buildPAN}}, \code{\link[igraph:layout]{layout}}
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
##view inferred PAN
viewPAN(pan, what='graph', layout="layout.fruchterman.reingold")	
}
}












