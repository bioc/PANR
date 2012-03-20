\name{data-Bakal2007}
\alias{Bakal2007}
\alias{Bakal2007Cluster}
\alias{nodeColor}
\alias{bm1}
\title{
Rich morphological phenotypes for gene overexpression and RNA 
interference screens
}                                                                        
\description{
The data set we use here comes from quantitative morphological 
screening for 249 gene-overexpression or RNAi knock-down 
experiments. For each individual cell, 145 different geometric 
features were computed by imaging analysis, and are subsequently 
scored with NNs trained to discriminate seven reference TCs 
with distinctive morphologies. For each TC, NN z-scores were 
computed from all scored cells in this TC (more details in 
Bakal 2007). 

\code{Bakal2007}: a matrix of NN z-scores with rows and columns 
corresponding to 249 TCs and seven reference TCs.

\code{Bakal2007Cluster}: unsupervised hierarchical clustering 
results by Bakal et al.

\code{nodeColor}: colors scaled according to the clustering 
results by Bakal et al.

\code{bm1}: an object of \code{\link[PAN:BetaMixture]{BetaMixture}}, 
which includes data and results for beta-mixture modelling 
on the association densities based on first-order cosine 
similarities of the phenotyping screens.
}
\usage{
##see example for details
}
\references{
Bakal, C. and Aach, J. and Church, G. and Perrimon, N. (2007). 
Quantitative morphological signatures define local signaling 
networks regulating cell morphology. Science, 316(5832), 1753.
}
\author{
Xin Wang \email{xw264@cam.ac.uk}
}

\examples{
data(Bakal2007)
dim(Bakal2007)
data(bm)
bm1
}










