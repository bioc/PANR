###########################1. scaling functions##############################
##A.scale node size
panNodeSizeScale<-function(x, minv=5, maxv=15) {
	minv+(maxv-minv)*(x-min(x))/(max(x)-min(x))
}
##B. scale node color
panNodeColScale<-function(x, n=10, neg.col="darkblue", pos.col="darkred", sep.col) {
	bneg<-0
	bpos<-0
	negColPl<-colorRampPalette(colors = c(neg.col, sep.col))
	posColPl<-colorRampPalette(colors = c(sep.col, pos.col))
	x.col<-rep("", length(x))
	names(x.col)<-names(x)
	##n means the number of color breaks in each half
	negCols<-negColPl(n+3)[1:(n+2)]
	x.col[x<bneg]<-negCols[as.integer(cut(x[x<bneg], breaks=c(min(x[x<bneg]-1), seq(min(x[x<bneg]), max(x[x<bneg]), 
											by=(max(x[x<bneg])-min(x[x<bneg]))/n), max(x[x<bneg])+1)))]
	posCols<-posColPl(n+3)[2:(n+3)]
	x.col[x>bpos]<-posCols[as.integer(cut(x[x>bpos], breaks=c(min(x[x>bpos]-1), seq(min(x[x>bpos]), max(x[x>bpos]), 
											by=(max(x[x>bpos])-min(x[x>bpos]))/n), max(x[x>bpos])+1)))]
	x.col[x==0]<-sep.col
	x.col
}
##C. scale edge width
panEdgeWidthScale<-function(x, minv=2, maxv=10) {
	minv+(maxv-minv)*(x-min(x))/(max(x)-min(x))
}
#########################2. build a RedeR graph from igraph###################
buildGraph4PAN<-function(pan, ig, engine, para) {
	#####################1. default settings for graph attributes################
	default.node.label.size<-ifelse(engine=="igraph", 10, 15)
	default.node.size.min<-ifelse(engine=="igraph", 5, 10)
	default.node.size.max<-ifelse(engine=="igraph", 20, 40)
	default.node.col<-gray(0.7)
#			default.node.low.col<-rgb(red=0, green=0, blue=139, maxColorValue=255)
#			default.node.high.col<-rgb(red=139, green=0, blue=0, maxColorValue=255)
	default.node.low.col<-colorRampPalette(c(colors()[90], colors()[90]))(1)
	default.node.high.col<-colorRampPalette(c(colors()[99], colors()[99]))(1)
	default.node.sep.col<-"#FFFFFF"
	default.edge.width.min<-1
	default.edge.width.max<-ifelse(engine=="igraph", 5, 15)
	default.edge.col<-gray(0.7)
#			default.edge.low.col<-rgb(red=0, green=139, blue=0, maxColorValue=255)
#			default.edge.high.col<-rgb(red=255, green=64, blue=64, maxColorValue=255)
	default.edge.low.col<-colorRampPalette(c("red", "red"))(1)
	default.edge.high.col<-colorRampPalette(c(colors()[51], colors()[51]))(1)
	default.gene.cluster.group.size<-5
	default.noOfModulePerRow<-4
	##pre-compute
	if(para$nodeSumMethod=="none")
		nodePheno<-pan@bm1@pheno[, para$nodeSumCols]
	else if(para$nodeSumMethod=="average")
		nodePheno<-rowSums(pan@bm1@pheno[, para$nodeSumCols])
	else if(para$nodeSumMethod=="median")
		nodePheno<-apply(pan@bm1@pheno[, para$nodeSumCols],1, median)
	else 
		stop("PAN-error: no other summarization method supported\n")
	#############################(1) Node##################################
	###A. node size
	ig.deg<-igraph::degree(ig, v=V(ig), mode = c("all"), loops = FALSE)			##compute degrees of nodes
	if(is.null(para$nodeSize)) 
		ig.nodeSize<-panNodeSizeScale(ig.deg, minv=default.node.size.min, maxv=default.node.size.max)
	else 
		ig.nodeSize<-para$nodeSize
	###B. node color
	if(is.null(para$nodeColor)) {
		ig.nodeColor<-default.node.col
		ig.nodeColor<-panNodeColScale(nodePheno[match(V(ig)$name, rownames(pan@bm1@pheno))], 
				neg.col=default.node.low.col, pos.col=default.node.high.col, sep.col=default.node.sep.col)
	} else {
		ig.nodeColor<-para$nodeColor
	}
	ig.nodeLineColor<-"#FFFFFF"
	ig.nodeFontSize<-default.node.label.size
	#############################(2) Edge##################################
	E(ig)$name<-paste("l", 1:length(E(ig)$weight),sep="")
	##A. edge width
	if(is.null(para$edgeWidth))
		ig.edgeWidth<-panEdgeWidthScale(E(ig)$weight, minv=default.edge.width.min, maxv=default.edge.width.max)
	else 
		ig.edgeWidth<-para$edgeWidth
	##B. edge color
	if(is.null(para$edgeColor)) {
		ig.edgeColor<-rep(default.edge.col, length(E(ig)$weight))
		ig.edgeColor[sign(E(ig)$weight)==1]<-default.edge.high.col
		ig.edgeColor[sign(E(ig)$weight)==-1]<-default.edge.low.col	
	} else {
		ig.edgeColor<-para$edgeColor
	}		
	##assign to igraph
	if(engine=="igraph") {
		V(ig)$size<-ig.nodeSize
		V(ig)$color<-ig.nodeColor
		V(ig)$frame.color<-ig.nodeLineColor
		V(ig)$label<-V(ig)$name
		E(ig)$width<-ig.edgeWidth
		E(ig)$color<-ig.edgeColor
	} else if(engine=="RedeR") {
		V(ig)$nodeSize<-ig.nodeSize
		V(ig)$nodeColor<-ig.nodeColor
		V(ig)$nodeLineColor<-ig.nodeLineColor
		V(ig)$nodeFontSize<-ig.nodeFontSize
		E(ig)$edgeWidth<-ig.edgeWidth
		E(ig)$edgeColor<-ig.edgeColor
	}
	#############################(3) Legend##################################
	legd<-list()
	##node size
	min2maxig.deg<-min(ig.deg):max(ig.deg)
	legd$nodeSize<-data.frame(degree=min2maxig.deg, diameter=panNodeSizeScale(min2maxig.deg, default.node.size.min, default.node.size.max))
	##node color
	min2maxig.zscore<-c(seq(min(nodePheno[nodePheno<0]), max(nodePheno[nodePheno<0]), by=(max(nodePheno[nodePheno<0])-min(nodePheno[nodePheno<0]))/10), 0, 
			seq(min(nodePheno[nodePheno>0]), max(nodePheno[nodePheno>0]), by=(max(nodePheno[nodePheno>0])-min(nodePheno[nodePheno>0]))/10))
	legd$nodeColor<-data.frame(zscore=min2maxig.zscore, color=panNodeColScale(min2maxig.zscore,neg.col=default.node.low.col, pos.col=default.node.high.col,sep.col="#FFFFFF"))
	##edge width
	min2maxig.SNR<-seq(min(E(ig)$weight), max(E(ig)$weight), by=(max(E(ig)$weight)-min(E(ig)$weight))/9)
	legd$edgeWidth<-data.frame(SNR=min2maxig.SNR, width=panEdgeWidthScale(min2maxig.SNR, minv=default.edge.width.min, maxv=default.edge.width.max))
	return(list(graph=ig, legend=legd))
}






































