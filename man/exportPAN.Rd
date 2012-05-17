\name{exportPAN}
\alias{exportPAN}
\alias{exportPAN,PAN,character_Or_missing,character_Or_missing,numeric_Or_integer_Or_missing,character_Or_missing,logical_Or_missing-method}
\title{
Export inferred PAN or module graphs to files
}
\description{
Powered by function \code{\link[igraph0:write.graph]{write.graph}} in package \code{\link[igraph0:igraph]{igraph}}, this function writes the 
inferred PAN or module graphs to files in a variety of formats that are supported 
by igraph0. 
}
\usage{
exportPAN(object, file="pan", what="graph", moduleID=1, format="gml", 
verbose=TRUE, ...)
}
\arguments{
	\item{object}{
an object of S4 class \code{PAN}.
}
	\item{file}{
the name of the file to write the graph or module (no need to put a filename 
extension).
}
	\item{what}{
a character value specifying what to write: `graph' or `pvclustModule'.
}
	\item{moduleID}{
a numeric or integer vector of modules to write (see details in \code{\link[PANR:sigModules]{sigModules}}).
}
	\item{format}{
a character value specifying the format to write (see more details in \code{\link[igraph0:write.graph]{write.graph}}).
}
	\item{verbose}{
a logical value to switch on (if \code{TRUE}) or off {if \code{FALSE}} detailed
run-time message.
}
	\item{...}{
not in use, but only for further extension.
}
}
\value{
a numeric vector of ids for significant gene modules
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
exportPAN(pan, file="pan", what="graph", format="ncol")
}
}












