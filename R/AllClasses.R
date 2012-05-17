###############################################################################
# All S4 classes for PAN
###############################################################################
##Class Beta-Mixture
setClass(
		"BetaMixture",
		representation(
				pheno="matrix",
				metric="character",
				order="numeric",
				association="numeric",
				model="character",
				partition="numeric",
				summary="list",
				result="list"
		),
		prototype=list(
				pheno=matrix(numeric(0), 0,0),
				metric="cosine",
				order=1,
				association=numeric(), 
				model="global",
				partition=integer(),
				summary=list(),
				result=list()
		)
)
setOldClass("igraph")
##Class PAN (Posterior association network)
setClass(
		"PAN",
		representation(
				bm1="BetaMixture",
				bm2="BetaMixture",
				edgeWt="matrix",
				graph="matrix",
				modules="list", 
				iPAN="igraph",
				legend="list",
				engine="character",
				summary="list"
		),
		prototype=list(
				bm1=new("BetaMixture"), 
				bm2=new("BetaMixture"), 
				edgeWt=matrix(),
				graph=matrix(),
				modules=list(),
				iPAN=igraph0::graph.empty(),
				legend=list(),
				engine="RedeR",
				summary=list()
		)
)

























