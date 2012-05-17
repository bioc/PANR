\name{viewLegend}
\alias{viewLegend}
\alias{viewLegend,PAN,character_Or_missing-method}
\title{
View the legends for the graph built for PAN
}
\description{
This function helps the user to view legends for the built PAN graph.
}
\usage{
viewLegend(object, what="nodeColor", ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{what}{
a character value specifying which legend to view: `nodeColor', `nodeSize'
or `edgeWidth'.
}
	\item{...}{
not in use, but only for further extension.
}
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
pan<-buildPAN(pan, engine="igraph0", para=list(nodeColor=nodeColor, 
hideNeg=TRUE), verbose=TRUE)
viewLegend(pan, "edgeWidth")
}
}













