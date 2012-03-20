##S4 methods for PAN
##initialization method
setMethod("initialize",
		c("PAN"),
		function(.Object, bm1, bm2) {
			#######################
			##check input arguments
			#######################
			if(is.null(bm1@result$fitBM))
				stop("PAN-error: No results found in 'bm1'")
			if(!missing(bm2) && is.null(bm2@result$fitBM))
				stop("PAN-error: No results found in 'bm2'")
			
			#######################
			##	initialization
			#######################
			##initialize log signal-noise ratios
#			object@iPAN<-igraph::graph.empty(n=0, directed=FALSE)
##			.Object@engine<-"igraph"
			if(!missing(bm2))
				.Object@bm2<-bm2
			.Object@bm1<-bm1
			nrR<-nrow(bm1@pheno)
			geneNames<-rownames(bm1@pheno)
			##
##			.Object@summary<-list()
			.Object@summary$input<-
					c(
							bm1=ifelse(missing(bm1), "None", "Y"),
							bm2=ifelse(missing(bm2), "None", "Y")
					)
			.Object@summary$graph<-NULL
			.Object@summary$pvclustModules<-NULL
			.Object@summary$PANModules<-NULL
			##return
			.Object
		}
)
##show summary information on screen
setMethod(
		"show",
		"PAN",
		function(object) {
			cat("A posterior association network:\n")
			summarize(object, what=c("input", "graph"))
		}
)
##print summary information on screen
setMethod(
		"summarize",
		c("PAN","character_Or_missing"),
		function(object, what="ALL") {
			##
			if(any(c("ALL","input") %in% what)) {
				cat("\n")
				cat("-input \n")
				if(!is.null(object@summary$input))
					print(object@summary$input, quote=FALSE, justify="centre")
				else
					cat("not available yet!")
				cat("\n")
			}
			if(any(c("ALL","graph") %in% what)) {
				cat("-graph \n")
				if(!is.null(object@summary$graph))
					print(object@summary$graph, quote=FALSE, justify="centre")
				else
					cat("not available yet!")
				cat("\n")
			}
			if(any(c("ALL","module") %in% what)) {
				cat("-modules \n")
				cat("--modules found by pvclust \n")
				if(!is.null(object@summary$pvclustModules))
					print(object@summary$pvclustModules, quote=FALSE, justify="centre")
				else
					cat("not available yet!\n")
			}
		}
)
##compute edge weights
##type: SNR, PPR, PP, 
setMethod(
		"edgeWeight",
		c("PAN", "character_Or_missing", "character_Or_missing", "logical_Or_missing"),
		function(object, which="bm1", type="SNR", log=TRUE, ...) {
			if(which=="bm1") {
				pZ<-object@bm1@result$fitBM$z
				nrGene<-nrow(object@bm1@pheno)
				geneNames<-rownames(object@bm1@pheno)
			} else if(which=="bm2") {
				pZ<-object@bm2@result$fitBM$z
				nrGene<-nrow(object@bm2@pheno)
				geneNames<-rownames(object@bm2@pheno)
			}			
			##snr
			if(type=="SNR") 
				wt<-((pZ[, 1]+pZ[, 3])/pZ[, 2])
			##posterior probability
			else if(type=="PPR") {
				wt<-pZ[, 3]/(1-pZ[, 3])	
			}
			else if(type=="PP") {
				wt<-pZ[, 3]
			}
			if(log)
				wt<-log(wt)
			##transform back to matrix
			wtMat<-matrix(0, nrGene, nrGene)
			dimnames(wtMat)<-list(geneNames, geneNames)
			wtMat[upper.tri(wtMat)]<-wt	
			wtMat<-wtMat+t(wtMat)			##make it symmetric
			wtMat
		}
)

##inference
setMethod(
		"infer",
		c("PAN", "list", "logical_Or_missing", "logical_Or_missing"),
		function(object, para=list(type="SNR", log=TRUE, sign=TRUE, cutoff=0), filter=FALSE, verbose=TRUE, ...) {		
			##1 check parameters
			def_paraNames<-c("type", "log", "sign", "cutoff")
			paraNames<-names(para)
			ukn_paraNames<-setdiff(paraNames, def_paraNames)
			if(length(ukn_paraNames)>0)
				warning("PAN-warning: input unknown parameters: ", paste(ukn_paraNames,collapse=", "))
			###check type
			if(is.null(para$type)) {
				para$type<-"SNR"
			} else if(!para$type%in%c("SNR", "PPR", "PP"))
				stop("PAN-error: para$type must be 'SNR', 'PPR' or 'PP'")
			###check log
			if(!is.null(para$log)) {
				if(!is.logical(para$log))
					stop("PAN-error: para$log must be a logical value")
			} else 
				para$log<-TRUE
			###check sign
			if(!is.null(para$sign)) {
				if(!is.logical(para$sign))
					stop("PAN-error: para$sign must be a logical value")
			} else 
				para$sign<-TRUE
			###check cutoff
			if(!is.null(para$cutoff)) {
				if(!is.numeric(para$cutoff) && !is.numeric(para$cutoff) && para$cutoff<0)
					stop("PAN-error: para$cutoff must be a positive numeric value")
			} else 
				para$cutoff<-0
			###check type and sign together
			if(para$type=="PP" && (para$log || para$sign))
				stop("PAN-error: when para$type='PP', para$log and para$sign must be both FALSE")
			if(para$type=="PPR" && (para$sign))
				stop("PAN-error: when para$type='PPR', para$sign must be FALSE")
			##2 compute edge weights
			wtMat<-edgeWeight(object, which="bm1", type=para$type, log=para$log)	
			##find the index of the most likely edge model for each edge
			nrAsso<-nrow(object@bm1@result$fitBM$z)
			nrGene<-nrow(object@bm1@pheno)
			inds<-sapply(1:nrAsso, function(i) {
				which.max(object@bm1@result$fitBM$z[i,])
			}) 
			##sign
			if(para$sign) {
				signMat<-matrix(0, nrGene, nrGene)
				signMat[upper.tri(signMat)]<-sapply(1:nrAsso, function(i) {
							c(-1, 0, 1)[inds[i]]
						})
				signMat<-signMat+t(signMat)
			}
			if(para$sign)
				object@edgeWt<-wtMat*signMat
			else 
				object@edgeWt<-wtMat
			##infer a PAN based on cutoff of edge weights
			object@graph<-object@edgeWt
			object@graph[which(wtMat<para$cutoff, arr.ind=TRUE)]<-0			
			if(sum(object@graph)==0)
				warning("PAN-warning: no edge is found significant based on the specified 'snrCutoff'")
			else if(filter) {
				inds<-which(colSums(abs(object@graph))==0 & rowSums(abs(object@graph))==0)
				if(length(inds)>0)
					object@graph<-object@graph[-inds, -inds]
				nf<-nrow(object@edgeWt)-nrow(object@graph)
				if(verbose)
					if(nf>0)
						cat("-",nf, " genes have no significant interactions with any other gene \n")
			}
			##update summary
			object@summary$graph<-c(edgeWeightType=para$type, log=para$log, signed=para$sign, 
					edgeWeightCutoff=para$cutoff, nodes=nrow(object@graph), edges=sum(object@graph[upper.tri(object@graph)]>0))
			object
		}
)
##search enriched gene modules
setMethod(
		"pvclustModule",
		c("PAN", "numeric_Or_integer_Or_missing", "character_Or_missing", "character_Or_missing", "logical_Or_missing", "logical_Or_missing"),
		function(object, nboot=1000, metric="cosine", hclustMethod="average", filter=TRUE, verbose=TRUE, ...) {
			if(verbose)
				cat("-Searching enriched functional gene modules using clustering with bootstrap resampling\n")
			##use phenotypes without genes that have no significant associations with any other genes
			if(is.null(object@graph))
				stop("PAN-error: please do PAN inference before searching modules\n")
			if(filter)
					fPheno<-object@bm1@pheno[match(rownames(object@graph), rownames(object@bm1@pheno)), ]
				else
					fPheno<-object@bm1@pheno
			if(metric=="cosine")
				metric<-"uncentered"
			if(is(getOption("cluster"), "cluster") && 
					"package:snow" %in% search()) {
#				##!pach pvclust
#				ns<-getNamespace("pvclust")
#				en<-as.environment("package:pvclust")
#				
#				unlockBinding("dist.pvclust", ns)
#				assignInNamespace("dist.pvclust",dist.pvclust4PAN,ns="pvclust", envir=ns)
#				lockBinding("dist.pvclust", ns)
#				dist.pvclust<-getFromNamespace("dist.pvclust", ns=getNamespace("pvclust"))
#				
#				##!patch parPvclust
#				unlockBinding("parPvclust", ns)
#				assignInNamespace("parPvclust",parPvclust4PAN,ns="pvclust", envir=ns)
#				lockBinding("parPvclust", ns)
#				
#				unlockBinding("parPvclust", en)
#				assignInNamespace("parPvclust",parPvclust4PAN,ns="pvclust", envir=en)
#				lockBinding("parPvclust", en)			
#				##parPvclust<-getFromNamespace("parPvclust", ns)

				pvclust.rslt<-pvclust::parPvclust(getOption("cluster"), t(fPheno), method.hclust=hclustMethod,method.dist=metric, nboot=nboot, ...)
			} else {
				pvclust.rslt<-pvclust(t(fPheno), method.hclust=hclustMethod,method.dist=metric, nboot=nboot, ...)
			}
			if(is.null(object@modules))
				object@modules<-list()
			object@modules$pvclust<-pvclust.rslt
			object@modules$pval<-1-pvclust.rslt$edges[, "au"]
			##
			hc2split4PAN<-getFromNamespace("hc2split", "pvclust")
			##
			object@modules$clusters<-lapply(hc2split4PAN(pvclust.rslt$hclust)$member, function(x) pvclust.rslt$hclust$label[x])
			##
			object@summary$pvclustModules<-c(nboot=nboot, hclustLinkage=hclustMethod, "significantModules(pval<0.05)"=sum(object@modules$pval<0.05))
			object
		}
)

##search enriched gene modules
setMethod(
		"exportPAN",
		c("PAN", "character_Or_missing", "character_Or_missing", "numeric_Or_integer_Or_missing", "character_Or_missing", "logical_Or_missing"),
		function(object, file="pan", what="graph", moduleID=1, format="gml", verbose=TRUE, ...) {
			if(is.null(object@iPAN))
				stop("PAN-error: please build igraph first using 'buildPANGraph'")
			########################0. check validity of arguments ######################
			##A. check 'what'
			if(!what%in%c("graph", "pvclustModule"))
				stop("PAN-error: 'what' should be either 'graph' or 'pvclustModule'")
			##B. check 'moduleID'
			if(what=="pvclustModule" &&(any(moduleID<0) || any(moduleID>length(object@modules$clusters))))
				stop(paste("PAN-error: 'moduleID' should be >0 and <", length(object@modules$clusters), sep=""))
			##C. check 'format'
			if(!format%in%c('pajek', 'graphml', 'dot', 'gml', 'edgelist', 'lgl', 'ncol', 'dimacs'))
				stop("PAN-error: 'format' should be: 'pajek', 'graphml', 'dot', 'gml', 'edgelist', 'lgl', 'ncol', or 'dimacs'")
			########################1. export a PAN or modules###########################
			if(verbose)
				cat(paste("-exporting ", ifelse(what=='graph', "a PAN", "PAN modules"), " ...\n",sep=""))
			if(what=="graph") {
				igraph::write.graph(object@iPAN, file=paste(file, ".", format,sep=""), format=format)
			} else if(what=="pvclustModule") {
				##reorder according to pval
				sigmodulesWithInRange<-object@modules$clusters[moduleID]
				sigmodulesPval<-object@modules$pval[moduleID]
				##write tables
				for(mod in 1:length(sigmodulesWithInRange)) {
					sub.ig<-subgraph(object@iPAN, sigmodulesWithInRange[[mod]])
					igraph::write.graph(sub.ig, file=paste(file, "_module", mod, "(",format(sigmodulesPval[mod], digits=3, nsmall=2, scientific=TRUE), ").", format,sep=""), format=format)
				}
			}
			if(verbose)
				cat(paste("-finished exporting ", ifelse(what=='graph', "a PAN", "PAN modules"), " ...\n",sep=""))
		}
)
##sigModules
setMethod(
		"sigModules",
		c("PAN", "numeric_Or_integer_Or_missing", "numeric_Or_integer_Or_missing", "numeric_Or_integer_Or_missing", "character_Or_missing", "logical_Or_missing"),
		function(object, pValCutoff=0.01, minSize=3, maxSize=100, sortby="size", decreasing=FALSE, ...) {
			if(is.null(object@modules))
				stop("PAN-error: please search modules first using 'pvclustModule'")
			if(pValCutoff<0 || minSize<0 || maxSize<0)
				stop("PAN-error: input arguments should be >0")
			if(!sortby%in%c("pval", "size"))
				stop("PAN-error: 'sortby' can only be 'pval' or 'size'")
			lclusSize<-unlist(lapply(object@modules$clusters, length))
			sigInds<-which(object@modules$pval<pValCutoff & lclusSize>=minSize & lclusSize<=maxSize)
			if(length(sigInds)>0) {
				if(sortby=="pval") {
					sigInds<-sigInds[sort.list(object@modules$pval[sigInds], decreasing=decreasing)]
				} else if(sortby=="size") {
					sigInds<-sigInds[sort.list(lclusSize[sigInds], decreasing=decreasing)]
				}
			}
			sigInds
		}
)
##view nested modules in a graph
setMethod(
		"viewNestedModules", 
		c("PAN", "numeric_Or_integer_Or_missing", "numeric_Or_integer_Or_missing", "numeric_Or_integer_Or_missing", "logical_Or_missing"),
		function(object, pValCutoff=0.01, minSize=3, maxSize=100, verbose=TRUE, ...) {
			if(is.null(object@iPAN))
				stop("PAN-error: please build igraph first using 'buildPANGraph'")
			if(object@engine!="RedeR")
				stop("PAN-error: the current version can only view the nested modules in 'RedeR'")
			##		
			sigInds<-sigModules(object, pValCutoff, minSize, maxSize)
			moduleNameList<-object@modules$clusters[sigInds]
			moduleNameList<-moduleNameList[sort.list(unlist(lapply(moduleNameList, length)), decreasing=TRUE)]
			if(length(sigInds)==0)
				stop("PAN-error: no significant modules found given the arguments ")
			
			moduleRef<-list()
			isParents<-list()
			globalMods<-c(1)			
			for(i in 1:length(moduleNameList)) {
				##find parents
				if(i>1) {
					isParents[[i]]<-sapply(1:(i-1), function(pt) {
								if(all(moduleNameList[[i]]%in%moduleNameList[[pt]]))
									return(TRUE)
								else
									return(FALSE)
							})
					if(sum(isParents[[i]])==0)
						globalMods<-c(globalMods, i)
				}
			}
			default.grid.width<-22
			default.grid.height<-30
			default.noOfModulesPerRow<-4
			nrGlobalMods<-length(globalMods)
			nrModRows<-ceiling(nrGlobalMods/default.noOfModulesPerRow)
			modCoordYs<-(1:nrModRows-1)*default.grid.height+15
			modCoordXs<-(1:default.noOfModulesPerRow-1)*default.grid.width+12.5
			tempMat<-matrix(c(globalMods, rep(NA, default.noOfModulesPerRow*nrModRows-nrGlobalMods)), nrow=nrModRows, ncol=default.noOfModulesPerRow, byrow=TRUE)
			
			nrNodes<-length(V(object@iPAN))
			##upload graph
			rdp <- RedPort('MyPort')
			calld(rdp)
			object@iPAN$zoom<-70
			gigRef<-addGraph(rdp, object@iPAN, layout.random(object@iPAN))
			for(i in 1:length(moduleNameList)) {
				ml<-length(moduleNameList[[i]])
				if(i>1 && sum(isParents[[i]])>0) {							
					closestParent<-rev(which(isParents[[i]]))[1]
					moduleRef[[i]]<-nestNodes(rdp, nodes=moduleNameList[[i]], parent=moduleRef[[closestParent]], 
							gscale=100*ml/length(moduleNameList[[closestParent]]), gatt=list(isAnchor=FALSE, nestLineType="DOTTED", nestShape="ROUNDED_RECTANGLE", 
									nestLineWidth=3, nestColor="#FFFFFF", nestLineColor="#000000"))
				} else {
					xi<-which(tempMat==i, arr.ind=TRUE)
					gscale<-round(100*ml/nrNodes)
					if(gscale<30)
						gscale<-30
					moduleRef[[i]]<-nestNodes(rdp, nodes=moduleNameList[[i]], parent=gigRef, 
							gcoord=c(modCoordXs[xi[1, 2]], modCoordYs[xi[1, 1]]), gatt=list(isAnchor=TRUE, nestLineType="DOTTED", nestShape="ROUNDED_RECTANGLE", 
									nestLineWidth=8, nestColor="#FFFFFF", nestLineColor="#000000"),gscale=gscale)
				}
			}
			mergeOutEdges(rdp, lb=1, ub=15)
		}
)
##view in igraph or RedeR
setMethod(
		"viewPAN",
		c("PAN", "character_Or_missing", "numeric_Or_integer_Or_missing", "character_Or_missing", "logical_Or_missing"),
		function(object, what="graph", moduleID=1, layout="layout.fruchterman.reingold", verbose=TRUE, ...) {
			if(is.null(object@iPAN))
				stop("PAN-error: please build igraph first using 'buildPANGraph'")
			########################0. check validity of arguments ######################
			##A. check 'what'
			if(!what%in%c("graph", "pvclustModule"))
				stop("PAN-error: 'what' should be either 'graph' or 'pvclustModule'")
			##B. check 'moduleID'
			if(what=="pvclustModule" &&(any(moduleID<0) || any(moduleID>length(object@modules$clusters))))
				stop(paste("PAN-error: 'moduleID' should be >0 and <", length(object@modules$clusters), sep=""))
			##C.
			if(object@engine=="igraph" && what=="pvclustModule" && length(moduleID)>1)
				warning("PAN-warning: multiple modules can only be viewed one after another in 'igraph' in the current version!")
			########################1. view graph in selected engine#####################
			pval.formatted<-format(object@modules$pval[moduleID], digits=3, nsmall=2, scientific=TRUE)
			if(object@engine=="RedeR") {
				rdp <- RedPort('MyPort')
				calld(rdp)
				if(what=="graph")
					addGraph(rdp, object@iPAN, layout=eval(parse(text=paste(layout,"(object@iPAN)",sep=""))))
				else if(what=="pvclustModule") {
					default.grid.width<-22
					default.grid.height<-30
					default.noOfModulesPerRow<-4
					nrGlobalMods<-length(moduleID)
					nrModRows<-ceiling(nrGlobalMods/default.noOfModulesPerRow)
					modCoordYs<-(1:nrModRows-1)*default.grid.height+15
					modCoordXs<-(1:default.noOfModulesPerRow-1)*default.grid.width+12.5
					tempMat<-matrix(c(moduleID, rep(NA, default.noOfModulesPerRow*nrModRows-nrGlobalMods)), nrow=nrModRows, ncol=default.noOfModulesPerRow, byrow=TRUE)
					for(i in 1:length(moduleID)) {
						xi<-which(tempMat==moduleID[i], arr.ind=TRUE)
						sub.ig<-subgraph(object@iPAN, object@modules$cluster[[moduleID[i]]])
						sub.ig$isNest<-TRUE
						sub.ig$isAssign<-TRUE
						sub.ig$coordX<-modCoordXs[xi[1, 2]]
						sub.ig$coordY<-modCoordYs[xi[1, 1]]
						sub.ig$isAnchor<-TRUE
						sub.ig$nestAliases<-paste("module", i, " (p.value=", pval.formatted[i],")",sep="")
						sub.ig$nestFontSize<-12
						sub.ig$nestSize<-250
						sub.ig$nestLineWidth<-3
						sub.ig$nestShape<-"ROUNDED_RECTANGLE"
						sub.ig$nestColor<-"#FFFFFF"
						sub.ig$nestLineColor<-"#000000"
						sub.ig$nestLineType<-"DOTTED"
						sub.ig$nestFontX<-9
						sub.ig$nestFontY<-12
						sub.ig$zoom<-50
						sub.ig$nestFontSize<-20
						sub.ig$nestLineWidth<-4
						addGraph(rdp, sub.ig, gscale=25, layout=eval(parse(text=paste(layout,"(sub.ig)",sep=""))))
					}
				}
			} else if(object@engine=="igraph") {
				if(what=="graph")
					plot.igraph(object@iPAN, layout=eval(parse(text=layout)))
				else if(what=="pvclustModule") {
					for(i in 1:length(moduleID)) {
						sub.ig<-subgraph(object@iPAN, object@modules$cluster[[moduleID[i]]])
						plot.igraph(sub.ig, layout=eval(parse(text=layout)), main=paste("module", i, " (p.value=", pval.formatted[i],")",sep=""))
						if(i<length(moduleID))
							readline(prompt = "(Press <Enter> to view the next module ...)") 
					}
				}
			}
		}
)
##build igraph or RedeR graph
setMethod(
		"buildPAN",
		c("PAN", "character_Or_missing", "list_Or_missing", "logical_Or_missing"),
		function(object, engine="igraph", para=list(nodeColor=NULL, nodeSize=NULL, edgeWidth=NULL, edgeColor=NULL, 
						nodeSumCols=1, nodeSumMethod="none", hideNeg=TRUE), verbose=TRUE, ...) {
			########################0. check validity of arguments ######################
			##A. check 'engine'
			if(!engine%in%c("igraph", "RedeR"))
				stop("PAN-error: 'engine' should be either 'igraph' or 'RedeR'")
			##B. check 'para'
			para0=list(nodeSumCols=1, nodeSumMethod="none", hideNeg=TRUE)
			moreparas<-setdiff(names(para0), names(para))
			if(length(moreparas)>0 && verbose)
				cat("--use default args: ", paste(moreparas, collapse=', '), "\n")
			para<-c(para0[moreparas], para)

			##graph
			if(is.null(object@graph))
				stop("PAN-error: please run 'infer' first")
			if(para$hideNeg) {
				posGraph<-object@graph
				posGraph[posGraph<0]<-0
				ig<-graph.adjacency(adjmatrix=posGraph, diag=FALSE, mode="undirected", weighted=TRUE)
			}
			else
				ig<-graph.adjacency(adjmatrix=(object@graph), diag=FALSE, mode="undirected", weighted=TRUE)
			if(verbose)
				cat(paste("-building ", engine, " graph ...\n",sep=""))
			ig<-buildGraph4PAN(object, ig, engine, para)
			if(verbose)
				cat(paste("-finished building ", engine, " graph ...\n",sep=""))
			object@iPAN<-ig$graph
			object@legend<-ig$legend
			object@engine<-engine
			object
		}
)
##build igraph or RedeR graph
setMethod(
		"viewLegend",
		c("PAN", "character_Or_missing"),
		function(object, what="nodeColor", ...) {
			if(is.null(object@legend))
				stop("PAN-error: please run `buildPAN' first to build a graph for inferred PAN ")
			if(!what%in%c("nodeColor", "edgeWidth", "nodeSize"))
				stop("PAN-error: 'what' should be: 'nodeColor', 'edgeWidth', 'nodeSize'")
			if(what=="edgeWidth") {
				##edge width
				edgesnr<-round(object@legend$edgeWidth[, "SNR"], 2)
				edgew<-object@legend$edgeWidth[, "width"]
				min.edgesnr<-min(edgesnr)
				max.edgesnr<-max(edgesnr)
				min.edgew<-min(edgew)
				max.edgew<-max(edgew)
				par(mar=c(6, 1, 6, 8))
				plot.new()
				plot.window(xlim=range(edgesnr), ylim=c(0, max.edgew))
				polygon(x=c(min.edgesnr, min.edgesnr, max.edgesnr, max.edgesnr), y=c(0, min.edgew, max.edgew, 0))
				axis(side=1, at=edgesnr, labels=format(edgesnr, ndigits=3, nsmall=2))
				mtext(side=1, line=3, text="log(signal-to-noise ratio)")
				axis(side=4, at=range(edgew), labels=round(range(edgew)))
				mtext(side=4, line=3, text="edge width")
			} else if(what=="nodeColor") {
				##node color
				plot.new()
				plot.window(xlim=c(1,nrow(object@legend$nodeColor)), ylim=c(0, 1))
				image(matrix(1:nrow(object@legend$nodeColor), ncol = 1), col = as.character(object@legend$nodeColor[, "color"]), axes=FALSE, xaxs="i",xlab="",ylab="",ylim=c(0,1))
				box()
				axis(side=1, at=seq(0,1,length=nrow(object@legend$nodeColor)), labels=format(round(object@legend$nodeColor[, 1], 2), nsmall=2))
				mtext(side=1, line=3, text="z-score")
				title(main="node color")
			} else if(what=="nodeSize") {
				##node size
				degs<-seq(min(object@legend$nodeSize[, "degree"]),max(object@legend$nodeSize[, "degree"]), by=4)
				degs<-unique(c(degs, max(object@legend$nodeSize[, "degree"])))
				nodesize<-object@legend$nodeSize[match(degs, object@legend$nodeSize[, "degree"]), ]
				legdIG<-matrix(0, nrow(nodesize), nrow(nodesize))
				rownames(legdIG)<-1:nrow(nodesize)
				colnames(legdIG)<-1:nrow(nodesize)
				graph.adjacency(adjmatrix=legdIG, )->legdIG
				V(legdIG)$size<-nodesize[, "diameter"]*0.4
				V(legdIG)$color<-"black"
				V(legdIG)$label<-""
				
				plot.igraph(legdIG, layout=cbind(seq(-1, 1, length=nrow(nodesize)), rep(-1, nrow(nodesize))))
				title(main="node size")
				axis(side=1, at=seq(-1, 1, length=nrow(nodesize)), labels=nodesize[, "degree"], tick=FALSE)
				mtext(side=1, line=3, text="degree")
			}
		}
)














