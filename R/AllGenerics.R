##generic functions
setGeneric("permNULL", 
	function(object, permMethod="keepRep", ...) 
		standardGeneric("permNULL"), 
package="PANR")
setGeneric("fitNULL", 
	function(object, nPerm=20, thetaNULL=c(alphaNULL=4, betaNULL=4), 
			sumMethod="median", permMethod="keepRep", verbose=TRUE, ...) 
		standardGeneric("fitNULL"), 
package="PANR")
setGeneric("fitBM", 
	function(object, para=list(zInit=NULL, 
			thetaInit=c(alphaNeg=2, betaNeg=4, alphaNULL=4, betaNULL=4, alphaPos=4, betaPos=2), gamma=NULL), 
			ctrl=list(fitNULL=FALSE, tol=1e-3, maxIter=NULL), verbose=TRUE, ...) 
		standardGeneric("fitBM"), 
package="PANR")
setGeneric("p2SNR", 
	function(object, pval, ...) 
		standardGeneric("p2SNR"), 
package="PANR")
setGeneric("SNR2p", 
	function(object, SNR, ...) 
		standardGeneric("SNR2p"), 
package="PANR")
setGeneric("infer", 
	function(object, para=list(type="SNR", log=TRUE, 
			sign=TRUE, cutoff=0), filter=FALSE, verbose=TRUE, ...) 
		standardGeneric("infer"), 
package="PANR")
setGeneric("edgeWeight", 
	function(object, which="bm1", type="SNR", log=TRUE, ...) 
		standardGeneric("edgeWeight"), 
package="PANR")
setGeneric("pvclustModule", 
	function(object, nboot=1000, metric="cosine", hclustMethod="average", 
		filter=TRUE, verbose=TRUE,  ...) 
		standardGeneric("pvclustModule"), 
package="PANR")
setGeneric("view", 
	function(object, what="fitBM", ...) 
		standardGeneric("view"), 
package="PANR")
setGeneric("exportPAN", 
		function(object, file="pan", what="graph", moduleID=1, format="gml", verbose=TRUE, ...)
			standardGeneric("exportPAN"), 
		package="PANR")
setGeneric("sigModules", 
		function(object, pValCutoff=0.01, minSize=3, maxSize=100, sortby="size", decreasing=FALSE, ...)
			standardGeneric("sigModules"), 
		package="PANR")
setGeneric("viewNestedModules", 
		function(object, pValCutoff=0.01, minSize=3, maxSize=100, verbose=TRUE, ...)
			standardGeneric("viewNestedModules"), 
		package="PANR")
setGeneric("viewPAN", 
		function(object, what="graph", moduleID=1, layout="layout.fruchterman.reingold", verbose=TRUE, ...)
			standardGeneric("viewPAN"), 
		package="PANR")
setGeneric("viewLegend", 
		function(object, what="nodeColor", ...)
			standardGeneric("viewLegend"), 
		package="PANR")
setGeneric("buildPAN", 
		function(object, engine="igraph", para=list(nodeColor=NULL, nodeSize=NULL, edgeWidth=NULL, edgeColor=NULL, 
						nodeSumCols=1, nodeSumMethod="none", hideNeg=TRUE), verbose=TRUE, ...)
			standardGeneric("buildPAN"), 
		package="PANR")
setGeneric("summarize",
	function(object, what="ALL", ...) 
		standardGeneric("summarize"), 
	package="PANR")
