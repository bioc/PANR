\name{PAN-class}
\alias{PAN}
\alias{PAN-class}
\docType{class}
\title{
An S4 class for inferring a posterior association network
}
\description{
This S4 class includes methods to infer posterior association networks
and enriched modules  of functional gene interactions from rich
phenotyping screens. 
}
\section{Objects from the Class}{ 
Objects of class \code{PAN} can be created from \code{new("PAN", bm1,
bm2)} (see the example below for details).
}
\section{Slots}{
	\describe{	
		\item{\code{bm1}:}{
		an object of S4 class \code{BetaMixture}, which models the first-
		order similarities between genes (see \code{\link[PANR:BetaMixture]{BetaMixture}}).
		}
		\item{\code{bm2}:}{
		an object of S4 class \code{BetaMixture}, which models the second-
		order similarities between genes (modularity).
		}
		\item{\code{edgeWt}:}{
		a weighted adjacency matrix computed from the posterior probabilities
		for gene associations to belong to mixture components (see
		\code{\link[PANR:edgeWeight]{edgeWeight}}).
		}
		\item{\code{engine}:}{
		the graphics visualization engine for PAN.
		}
		\item{\code{graph}:}{
		a weighted adjacency matrix with edge weights satisfying certain
		constraints specified by the user (see \code{\link[PANR:infer]{infer}}).
		}
		\item{\code{modules}:}{
		a list summarizing inferred enriched functional gene modules (see
		\code{\link[PANR:pvclustModule]{pvclustModule}}.
		}
		\item{\code{iPAN}:}{
		an igraph0 object for storing the inferred PAN.
		}
		\item{\code{legend}:}{
		a list of legends for built PAN graph.
		}
		\item{\code{summary}:}{
		a list of summary information for available results.
		}
	}
}
\section{Methods}{
	An overview of methods (More detailed introduction can be found in
	help for each specific function.):
	\describe{	
		\item{\code{edgeWeight}}{
		compute edge weights by signal-to-noise ratio, posterior odd or
		posterior probabilities (more details in \code{\link[PANR:edgeWeight]{edgeWeight}}).
		}
		\item{\code{infer}}{
		infer a posterior association network given the beta-mixture
		model(s) fitted to first- and/or second-order similarities (more
		details in \code{\link[PANR:infer]{infer}}). 
		}
		\item{\code{pvclustModule}}{
		search significantly enriched functional gene modules by hierarchical
		clustering with bootstrap resampling based on the package \code{pvclust}
		(more details in \code{\link[PANR:pvclustModule]{pvclustModule}}).
		}
		\item{\code{exportPAN}}{
		export the inferred PAN or modules to file(s) in a variety of 
		formats (more details in \code{\link[PANR:exportPAN]{exportPAN}}).
		}
		\item{\code{sigModules}}{
		retrieve significant gene modules that satisfy the given p-value 
		cutoff and module size range  (more details in \code{\link[PANR:sigModules]{sigModules}}).
		}
		\item{\code{viewNestedModules}}{
		view a nested structure for gene modules searched by hierarchical 
		clustering (more details in \code{\link[PANR:viewNestedModules]{viewNestedModules}}).
		}
		\item{\code{viewPAN}}{
		view the inferred PAN or modules in \code{\link[igraph0:igraph]{igraph}} or \code{\link[RedeR:RedeR]{RedeR}} (more 
		details in \code{\link[PANR:viewPAN]{viewPAN}}).
		}
		\item{\code{buildPAN}}{
		build a PAN graph for visualization in \code{\link[igraph0:igraph]{igraph}} or \code{\link[RedeR:RedeR]{RedeR}}
		(more details in \code{\link[PANR:viewPAN]{viewPAN}}).
		}
		\item{\code{viewLegend}}{
		View the legends for the graph built for PAN.
		}
		\item{\code{summarize}}{
		summarize results including input data and parameters, inferred
		graph and modules.
		}
	}
}
\author{Xin Wang \email{xw264@cam.ac.uk}}
\references{
Xin Wang, Mauro Castro, Klaas W. Mulder and Florian Markowetz, Posterior 
association networks and enriched functional gene modules inferred from 
rich phenotypic perturbation screens, in preparation.
}
\seealso{
  \code{\link[PANR:edgeWeight]{edgeWeight}}
  \code{\link[PANR:infer]{infer}}
  \code{\link[PANR:pvclustModule]{pvclustModule}}
  \code{\link[PANR:exportPAN]{exportPAN}}
  \code{\link[PANR:sigModules]{sigModules}}
  \code{\link[PANR:viewPAN]{viewPAN}}
  \code{\link[PANR:viewNestedModules]{viewNestedModules}}
  \code{\link[PANR:summarize]{summarize}}
}
\examples{
\dontrun{
data(bm, package="PANR")
##create an object of `PAN'
pan<-new("PAN", bm1=bm1)
##infer a PAN
pan<-infer(pan, para=list(type="SNR", log=TRUE, sign=TRUE, cutoff=log(5)),
	filter=FALSE, verbose=TRUE)
##build a PAN graph for RedeR, hide negative edges
##using colors scaled based on the clustering results from Bakal et al. 2007 
data(Bakal2007Cluster)
pan<-buildPAN(pan, engine="RedeR", para=list(nodeColor=nodeColor, hideNeg=TRUE))
##view PAN in RedeR
library(RedeR)
viewPAN(pan, what="graph")
##print a summary of results
summarize(pan, "ALL")
}
}
\keyword{classes}
