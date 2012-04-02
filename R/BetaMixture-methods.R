###############################################################################
# S4 methods for Beta-Mixture models
###############################################################################
##initialization method
setMethod("initialize",
		"BetaMixture",
		function(.Object, metric="cosine", order=1, association, model="global", pheno, partition) {
			#######################
			##check input arguments
			#######################
			if(!metric%in%c("cosine", "correlation", "user-defined"))
				stop("PAN-error: 'metric' should be 'cosine', 'correlation', 'user-defined'\n")
			if(!order%in%c(1,2))
				stop("PAN-error: 'order' should be 1 or 2 \n")
			if(!missing(association) && any(association>1 || association<0))#
				stop("PAN-error: 'association' should be within [0,1]\n")
			if(!missing(partition) && !missing(pheno))
				if(length(partition)!=choose(nrow(pheno),2))
					stop("PAN-error: the length of 'partition' should be the same as choose(nrow(pheno),2)\n")
			#######################
			##	initialization
			#######################
##			cat("-Initialization\n")
			if(!missing(pheno))
				.Object@pheno<-pheno
			.Object@metric<-metric
			.Object@model<-model
			.Object@order<-order
			if(!missing(partition)) 
				.Object@partition<-partition
			if(missing(association)) {
				if(!missing(pheno)) {
					cat("--Computing association scores ...\n")	
					if(order==1)
						.Object@association<-assoScore(pheno=t(pheno), metric=metric, verbose=FALSE)
					else if(order==2)
						.Object@association<-assoScore(assoScore(pheno=t(pheno), metric=metric, upperTri=FALSE, transform=FALSE, verbose=FALSE), metric=metric, verbose=FALSE)
				}
			}
			else 
				.Object@association<-association
			##summarys
##			.Object@summary<-list()
			.Object@summary$input<-
				c(
					phenotype=ifelse(missing(pheno), "None", paste(nrow(pheno), " x ", ncol(pheno), sep="")),
					partition=ifelse(missing(partition), "None", length(unique(partition))), 
					model=model,
					metric=metric
				)
			.Object@summary$fitNULL<-NULL
			.Object@summary$fitBM<-NULL
			##results
##			.Object@result<-list()
##			cat("-Initialization finished\n")
			.Object
		}
)
##NULL permutation
setMethod(
		"permNULL",
		c("BetaMixture", "character_Or_missing"),
		function(object, permMethod="keepRep", ...) {
			conds<-unique(colnames(object@pheno))
			nrCond<-length(conds)
			nrCol<-ncol(object@pheno)
			nrRep<-nrCol/nrCond
			nrRow<-nrow(object@pheno)
			##method means to keep positions of replicates or permute all data
			if(permMethod=="all") {
				perm.pheno<-matrix(object@pheno[sample(1:(nrRow*ncol(object@pheno)), nrRow*ncol(object@pheno))], nrRow, ncol(object@pheno))
			} else if(permMethod=="col") {
				perm.pheno<-sapply(1:nrCol, function(icol) object@pheno[sample(1:nrRow, nrRow), icol])
			} else if(permMethod=="row") {
				perm.pheno<-t(sapply(1:nrRow, function(irow) object@pheno[irow, sample(1:nrCol, nrCol)]))
				##rownames(perm.pheno)<-rownames(object@pheno)
			} else if(permMethod=="keepRep") {
				perm.ind<-matrix(sort.list(rep(sample(1:(nrCond*nrRow), (nrCond*nrRow)),each=nrRep)),nrCond*nrRep, nrRow)
				perm.pheno<-t(matrix(t(object@pheno)[perm.ind], nrCond*nrRep, nrRow))
			} else if(permMethod=="keepRepCol") {
				perm.pheno<-matrix(0, nrRow, nrCol)
				sapply(1:nrCond, function(cond) perm.pheno[, which(colnames(object@pheno)==conds[cond])]<<-object@pheno[sample(1:nrRow, nrRow), which(colnames(object@pheno)==conds[cond])])
				##dimnames(perm.pheno)<-dimnames(object@pheno)
			} 
			rownames(perm.pheno)<-rownames(object@pheno)
			colnames(perm.pheno)<-colnames(object@pheno)
			perm.pheno
		}
)
##fit the NULL distribution to permuted data
setMethod(
		"fitNULL",
		c("BetaMixture", "numeric_Or_integer_Or_missing", "numeric_Or_integer_Or_missing", "character_Or_missing", "character_Or_missing", "logical_Or_missing"),
		function(object, nPerm=20, thetaNULL=c(alphaNULL=4, betaNULL=4), sumMethod="median", permMethod="keepRep", verbose=TRUE, ...) {
			cat("-NULL fitting\n")
			alphaNULLs<-rep(0, nPerm)
			betaNULLs<-rep(0, nPerm)
			
			object@result$fitNULL<-list()
			object@result$fitNULL$perm.pheno<-list()
			##object@result$fitNULL$perm.asso<-list()
			object@result$fitNULL$fit<-list()
			if(verbose) {
				cat("--Permutation and fitting the central components\n")
				pb <- txtProgressBar(style=3)
			}
				
			for(i in 1:nPerm) {
				##permutations
				perm.pheno<-permNULL(object, permMethod)
				##compute cosine similarities for permuted data
				##perm.asso<-assoScore(t(perm.pheno), metric=object@metric, upperTri=TRUE, verbose=FALSE, ...)
				if(object@order==1)
					perm.asso<-assoScore(pheno=t(perm.pheno), metric=object@metric, upperTri=TRUE, verbose=FALSE)
				else if(object@order==2)
					perm.asso<-assoScore(assoScore(pheno=t(perm.pheno), metric=object@metric, upperTri=FALSE, transform=FALSE, verbose=FALSE), metric=object@metric, verbose=FALSE)
				##fit
				fitdat<-suppressWarnings(fitdistr(perm.asso, "beta", start=list(shape1=thetaNULL["alphaNULL"], shape2=thetaNULL["betaNULL"])))
				alphaNULLs[i]<-fitdat$estimate[1]
				betaNULLs[i]<-fitdat$estimate[2]
				##save permutation data
				object@result$fitNULL$perm.pheno[[i]]<-perm.pheno
				##object@result$fitNULL$perm.asso[[i]]<-perm.asso
				object@result$fitNULL$fit[[i]]<-fitdat
				if(verbose) 
					setTxtProgressBar(pb, i/nPerm)
			}
			if(verbose) {
				cat("\n-Permutation and fitting finished\n")
				close(pb)
			}

			##summarize and return
			if(sumMethod=="median") {
				object@result$fitNULL$thetaNULL<-c(alphaNULL=median(alphaNULLs), betaNULL=median(betaNULLs))
			} else if(sumMethod=="mean") {
				object@result$fitNULL$thetaNULL<-c(alphaNULL=mean(alphaNULLs), betaNULL=mean(betaNULLs))
			} else 
				stop("PAN-error: no other summarization methods supported for current version")
			##update summary
			object@summary$fitNULL<-c(
					"No. of perm"=nPerm,
					"summarize method"=sumMethod,
					"permutation method"=permMethod,
					"shape1 (start->fitted)"=paste(format(thetaNULL["alphaNULL"], digits=3), format(object@result$fitNULL$thetaNULL["alphaNULL"], nsmall=3, digits=3), sep="->"),
					"shape2 (start->fitted)"=paste(format(thetaNULL["betaNULL"], digits=3), format(object@result$fitNULL$thetaNULL["betaNULL"], nsmall=3, digits=3), sep="->")
			)
			object
		}
)
##maximum likelihood estimation or maximum a posterior
setMethod(
		"fitBM",
		c("BetaMixture", "list_Or_missing", "list_Or_missing", "logical_Or_missing"),
		function(object, para=list(zInit=NULL, thetaInit=c(alphaNeg=2, betaNeg=4, alphaNULL=4, betaNULL=4, alphaPos=4, betaPos=2), gamma=NULL), 
				ctrl=list(fitNULL=FALSE, tol=1e-3, maxIter=NULL), verbose=TRUE, ...) {
			##
			if(verbose) {
				cat("-Beta-Mixture Modelling\n")
				cat("(")
				cat(paste(ifelse(is.null(para$gamma), "ML", "MAP"), " estimation for ", object@model, " Beta Mixture Model based on EM algorithm) \n",sep=""))	
			}		
			############################initialization and checking input#########################
			##1 initialize partitions 
			N<-length(object@association)
			if(object@model=="global") 
				object@partition<-rep(1, N)
			parts<-unique(object@partition)
			partCounts<-sapply(parts, function(part) sum(object@partition==part))
			nrPart<-length(parts)
			##2 check parameters
			def_paraNames<-c("zInit", "thetaInit", "gamma")
			paraNames<-names(para)
			ukn_paraNames<-setdiff(paraNames, def_paraNames)
			if(length(ukn_paraNames)>0)
				warning("PAN-warning: input unknown parameters: ", paste(ukn_paraNames,collapse=", "))
			##3 check and initialize other parameters
			####Z
			if(!is.null(para$zInit)) {
				if(!is.matrix(para$zInit) || nrow(para$zInit)!=length(object@association) || ncol(para$zInit)!=3)
					stop("PAN-error: 'para$zInit' must be a N(No. of associations) x 3(No. of mixture components) matrix")
				Z<-para$zInit
				pZ<-Z
			} else {
				Z<-matrix(0, N, 3)
				pZ<-matrix(0, N, 3)	
				Z[which(object@association>=0 & object@association<0.33), 1]<-1
				Z[which(object@association>=0.33 & object@association<0.67), 2]<-1
				Z[which(object@association>=0.67 & object@association<=1), 3]<-1
				pZ<-Z
			}
			####theta
			if(!is.null(para$thetaInit)) {
				if(!(is.numeric(para$thetaInit) || is.integer(para$thetaInit)) || any(para$thetaInit<0))
					stop("PAN-error: 'para$thetaInit' must be a a numeric/integer vector of positive values")
				thetaNames<-c("alphaNeg", "betaNeg", "alphaNULL", "betaNULL", "alphaPos", "betaPos")
				if(length(para$thetaInit)!=6 || 
						!all(names(para$thetaInit)%in%thetaNames))
					stop("PAN-error: 'para$thetaInit' must include ", paste(thetaNames, collapse=', '))
				theta<-para$thetaInit
			} else 
				theta<-c(alphaPos=2, betaPos=4, alphaNULL=4, betaNULL=4, alphaNeg=4, betaNeg=2)
			####gamma	
			if(!is.null(para$gamma)) {
				if((nrow(para$gamma)!=nrPart || ncol(para$gamma)!=3))
					stop("PAN-error: 'para$gamma' should be a N(No. of partitions) X 3(No. of mixture components) matrix")
				gamma<-para$gamma
			} else 
				gamma<-matrix(1, nrPart, 3)
			####PI
			PI<-matrix(0, nrPart, 3)	
			PIInds<-sapply(object@partition, match, parts)	
			##4 check control arguments
			def_ctrlNames<-c("fitNULL", "tol", "maxIter")
			ctrlNames<-names(ctrl)
			ukn_ctrlNames<-setdiff(ctrlNames, def_ctrlNames)
			if(length(ukn_ctrlNames)>0)
				warning("PAN-warning: input unknown control arguments: ", paste(ukn_ctrlNames,collapse=", "))
			####fitNULL
			if(!is.null(ctrl$fitNULL)) {
				if(!(is.logical(ctrl$fitNULL) && length(ctrl$fitNULL)==1))
					stop("PAN-error: 'ctrl$fitNULL' must be a logical value")
				fitNULL<-ctrl$fitNULL
			} else 
				fitNULL<-FALSE
			####tol
			if(!is.null(ctrl$tol)) {
				if(!(is.numeric(ctrl$tol) || is.numeric(ctrl$tol)))
					stop("PAN-error: 'ctrl$tol' must be a numeric/integer value ")
				tol<-ctrl$tol
			} else 
				tol<-1e-3
			####maxIter	
			if(!is.null(ctrl$maxIter)) {
				if(!(is.numeric(ctrl$maxIter) || is.integer(ctrl$maxIter)))
					stop("PAN-error: 'ctrl$maxIter' must be a numeric/integer value")
				maxIter<-ctrl$maxIter
			} else 
				maxIter<-Inf
			
			
			##5 EM
			################################EM starts here#################################
			##nlm fun11315,ction for M-step, nullp as fitted NULL parameters
			BM_M_nlm<-function(p, y, Z, PI, nullp) {
				sum(sapply(1:length(parts), function(part) {
						-sum(
								Z * t(log(PI[part, ]) + 
										t(log(cbind(
													dbeta(y, shape1=p[1], shape2=p[2]), 
													dbeta(y, shape1=nullp[1], shape2=nullp[2]), 
													dbeta(y, shape1=p[3], shape2=p[4])
													)
											  )
										)
								)
						)				
					})) + (-sum((gamma-1)*log(PI)))
			}
			##nlm function for M-step, NULL not fitted
			BM_M_nlm2<-function(p, y, Z, PI) {
				sum(sapply(1:length(parts), function(part) {
						-sum(
								Z * t(log(PI[part, ]) + 
										t(log(cbind(
													dbeta(y, shape1=p[1], shape2=p[2]), 
													dbeta(y, shape1=p[3], shape2=p[4]), 
													dbeta(y, shape1=p[5], shape2=p[6])
													)
											  )
										)
								)
						)				
					})) + (-sum((gamma-1)*log(PI)))
			}
			NLL<-Inf
			irep<-0
			repeat{

				################################M-step#################################
				sapply(1:nrPart, function(i) {
					inds<-which(object@partition==parts[i])
					PI[i, ]<<-(colSums(pZ[inds, ])+gamma[i, ]-1)/(length(inds)+sum(gamma[i, ]-1))
				})
				if(fitNULL) {
					temp<-suppressWarnings(nlm(BM_M_nlm2, p=theta, y=object@association, Z=pZ, PI=PI, ...))
					theta<-temp$estimate
				}	else {
					temp<-suppressWarnings(nlm(BM_M_nlm, p=theta[c(1:2, 5:6)], 
							y=object@association, Z=pZ, PI=PI, nullp=theta[3:4], ...))
					theta<-c(temp$estimate[1:2], theta[3:4], temp$estimate[3:4])
				}
				################################E-step#################################	
				fmat<-PI[PIInds,]*(cbind(dbeta(object@association, shape1=theta[1], shape2=theta[2]), 
							dbeta(object@association, shape1=theta[3], shape2=theta[4]), 
							dbeta(object@association, shape1=theta[5], shape2=theta[6])))
				pZ<-fmat/apply(fmat, 1, sum)
				##stop
				if(abs(NLL-temp$minimum)<=tol || irep==maxIter) break
				NLL<-temp$minimum
				irep<-irep+1
				if(verbose)
					cat("--iter ", irep, "\t NLL ", round(NLL, 3), "\r")
			}	
			if(verbose) 
				cat("\n-Beta-Mixture Modelling fitting finished!\n")
			##6 update object and return
			##update result
			colnames(pZ)<-c("-", "x", "+")
			object@result$fitBM<-list()
			object@result$fitBM$z<-pZ
			object@result$fitBM$theta<-theta
			names(object@result$fitBM$theta)<-c("alphaNeg", "betaNeg", "alphaNULL", "betaNULL", "alphaPos", "betaPos")
			object@result$fitBM$pi<-PI
			object@result$fitBM$NLL<-NLL
			##update summary
			object@summary$fitBM<-list(
					"parameters"=
							c(
								"z (init)"=ifelse(is.null(para$zInit), "None", paste(dim(para$zInit), collapse=' x ')), 
								"shapes(- init)"=paste("shape1=",format(para$thetaInit["alphaNeg"], digits=3), 
														", shape2=", format(para$thetaInit["betaNeg"], digits=3), sep=""),
								"shapes(x init)"=paste("shape1=",format(para$thetaInit["alphaNULL"], digits=3), 
														", shape2=", format(para$thetaInit["alphaNULL"], nsmall=3, digits=3), sep=""),
								"shapes(+ init)"=paste("shape1=",format(para$thetaInit["alphaPos"], digits=3), 
														", shape2=", format(para$thetaInit["betaPos"], digits=3), sep=""),
								"gamma"=ifelse(is.null(para$gamma), "None", paste(dim(para$gamma), collpase=' x '))
							),
					"control args"=c(
								"fitNULL"=as.character(fitNULL),
								"tolerance"=as.character(tol),
								"maxIteration"=as.character(maxIter)
							),
					"result"=c(
							"shapes(- fitted)"=paste("shape1=",format(theta[1], nsmall=3, digits=3), 
														", shape2=", format(theta[2], nsmall=3, digits=3), sep=""),
							"shapes(x fitted)"=paste("shape1=",format(theta[3], nsmall=3, digits=3), 
														", shape2=", format(theta[4], nsmall=3, digits=3), sep=""),
							"shapes(+ fitted)"=paste("shape1=",format(theta[5], nsmall=3, digits=3), 
														", shape2=", format(theta[6], nsmall=3, digits=3), sep=""),
							"pi (fitted)"=paste(apply(format(PI, nsmall=3, digits=3), 1, paste, collapse=','), collapse=';'),
							"NLL"=as.character(format(NLL, nsmall=3, digits=3))
							)
					)
			object
		}
)
##translate p-value to SNR
setMethod(
		"p2SNR",
		c("BetaMixture", "numeric_Or_integer"),
		function(object, pval, ...) {
			theta<-object@result$fitBM$theta
			pis<-object@result$fitBM$pi
			x.l<-qbeta(p=pval, shape1=theta[3], shape2=theta[4], 
				lower.tail = TRUE, log.p = FALSE)
			x.h<-qbeta(p=pval, shape1=theta[3], shape2=theta[4], 
				lower.tail = FALSE, log.p = FALSE)
			if(object@model=="global") {
				snr.l<-(pis[1]*dbeta(x.l, shape1=theta[1], shape2=theta[2])+
					pis[3]*dbeta(x.l, shape1=theta[5], shape2=theta[6]))/
					(pis[2]*dbeta(x.l, shape1=theta[3], shape2=theta[4]))
				snr.h<-(pis[1]*dbeta(x.h, shape1=theta[1], shape2=theta[2])+
					pis[3]*dbeta(x.h, shape1=theta[5], shape2=theta[6]))/
					(pis[2]*dbeta(x.h, shape1=theta[3], shape2=theta[4]))
				rslt<-data.frame(pval=pval, q_lower=x.l, q_upper=x.h, 
					SNR_lower=snr.l, SNR_upper=snr.h)
			} else if(object@model=="global") {
				snr.l<-(pis[, 1]*dbeta(x.l, shape1=theta[1], shape2=theta[2])+
					pis[, 3]*dbeta(x.l, shape1=theta[5], shape2=theta[6]))/
					(pis[, 2]*dbeta(x.l, shape1=theta[3], shape2=theta[4]))
				snr.h<-(pis[, 1]*dbeta(x.h, shape1=theta[1], shape2=theta[2])+
					pis[, 3]*dbeta(x.h, shape1=theta[5], shape2=theta[6]))/
					(pis[, 2]*dbeta(x.h, shape1=theta[3], shape2=theta[4]))		
				rslt<-data.frame(pval=rep(pval, 2), q_lower=rep(x.l, 2), 
					q_upper=rep(x.h, 2), SNR_lower=snr.l, SNR_upper=snr.h)		
			}
			print(rslt)
			invisible(rslt)
		}
)
##translate SNR to p-value
setMethod(
		"SNR2p",
		c("BetaMixture", "numeric_Or_integer"),
		function(object, SNR, ...) {
			theta<-object@result$fitBM$theta
			pis<-object@result$fitBM$pi
			dbm<-function(x, lb, ub) {
				if(x<lb || x>ub) return(Inf)
				else 
				(pis[, 1]*dbeta(x, shape1=theta[1], shape2=theta[2])+
					pis[, 3]*dbeta(x, shape1=theta[5], shape2=theta[6]))/
					(pis[, 2]*dbeta(x, shape1=theta[3], shape2=theta[4]))
			}
			dbm2<-function(x, lb, ub, g=1) {
				if(x<lb || x>ub) return(Inf)
				else 
				(pis[g, 1]*dbeta(x, shape1=theta[1], shape2=theta[2])+
					pis[g, 3]*dbeta(x, shape1=theta[5], shape2=theta[6]))/
					(pis[g, 2]*dbeta(x, shape1=theta[3], shape2=theta[4]))
			}
			if(object@model=="global") {
				q.l<-suppressWarnings(nlm(f=function(x) {
					abs(SNR-dbm(x, 0, 0.5))}, p=0.01)$estimate)
				q.h<-suppressWarnings(nlm(f=function(x) {
					abs(SNR-dbm(x, 0.5, 1))}, p=0.99)$estimate)
				p.l<-pbeta(q=q.l, shape1=theta[3], shape2=theta[4], lower.tail = TRUE, 
					log.p = FALSE)	
				p.h<-pbeta(q=q.h, shape1=theta[3], shape2=theta[4], lower.tail = FALSE, 
					log.p = FALSE)	
				rslt<-data.frame(SNR=SNR, q_lower=q.l, q_upper=q.h, p_lower=p.l, 
					p_upper=p.h)
			} else if(object@model=="global") {
				q.l<-c(suppressWarnings(nlm(f=function(x) {
					abs(SNR-dbm2(x, 0, 0.5, 1))}, p=0.01)$estimate), 
					suppressWarnings(nlm(f=function(x) {
					abs(SNR-dbm2(x, 0, 0.5, 2))}, p=0.01)$estimate)
					)
				q.h<-c(suppressWarnings(nlm(f=function(x) {
					abs(SNR-dbm2(x, 0.5, 1, 1))}, p=0.99)$estimate), 
					suppressWarnings(nlm(f=function(x) {
					abs(SNR-dbm2(x, 0.5, 1, 2))}, p=0.99)$estimate)
					)
				p.l<-pbeta(q=q.l, shape1=theta[3], shape2=theta[4], lower.tail = TRUE, 
					log.p = FALSE)	
				p.h<-pbeta(q=q.h, shape1=theta[3], shape2=theta[4], lower.tail = FALSE, 
					log.p = FALSE)			
				rslt<-data.frame(SNR=SNR, q_lower=q.l, q_upper=q.h, p_lower=p.l, 
					p_upper=p.h)	
			}
			print(rslt)
			invisible(rslt)
		}
)
##BM density function
#setMethod(
#		"pdensity",
#		c("BetaMixture", "numeric_Or_integer"),
#		function(object, asso, ...) {
#			sum(PI*c(dbeta(x, shape1=p["alphaPos"], shape2=p["betaPos"]), dbeta(x, shape1=p["alphaNULL"], shape2=p["betaNULL"]), 
#							dbeta(x, shape1=p["alphaNeg"], shape2=p["betaNeg"])))
#		}
#)
##view fitting results
setMethod(
		"view",
		c("BetaMixture", "character_Or_missing"),
		function(object, what="fitBM", ...) {
			if(!(what %in% c("fitNULL", "fitBM"))) 
				step("PAN-error: 'what' should be 'fitNULL' or 'fitBM'")
			##model density
			betaMixture_overall_pdf<-function(PI, p, x) {
				sum(PI*c(dbeta(x, shape1=p["alphaNeg"], shape2=p["betaNeg"]), dbeta(x, shape1=p["alphaNULL"], shape2=p["betaNULL"]), 
								dbeta(x, shape1=p["alphaPos"], shape2=p["betaPos"])))
			}
			betaMixture_partial_pdf<-function(pi, a, b, x) {
				pi*dbeta(x, shape1=a, shape2=b)
			}
			##
			if(what=="fitNULL") {
				if(is.null(object@result$fitNULL))
					stop("PAN-error: 'fitNULL' result not available yet!")
				if(object@order==1)
					perm.asso<-assoScore(pheno=t(object@result$fitNULL$perm.pheno[[1]]), metric=object@metric, upperTri=TRUE, verbose=FALSE)
				else if(object@order==2)
					perm.asso<-assoScore(assoScore(pheno=t(object@result$fitNULL$perm.pheno[[1]]), metric=object@metric, upperTri=FALSE, 
									transform=FALSE, verbose=FALSE), metric=object@metric, verbose=FALSE)
				plot(density(perm.asso), col=gray(0.3), lwd=1, xlim=c(0,1), main="Fitting NULL distribution", xlab="Transformed cosine similarity")
				for(i in 2:length(object@result$fitNULL$perm.pheno)) {
					if(object@order==1)
						perm.asso<-assoScore(pheno=t(object@result$fitNULL$perm.pheno[[i]]), metric=object@metric, upperTri=TRUE, verbose=FALSE)
					else if(object@order==2)
						perm.asso<-assoScore(assoScore(pheno=t(object@result$fitNULL$perm.pheno[[i]]), metric=object@metric, upperTri=FALSE, 
										transform=FALSE, verbose=FALSE), metric=object@metric, verbose=FALSE)
					lines(density(perm.asso), col=gray(0.3), lwd=1, xlim=c(0,1))	
				}
				lines(1:49*0.02, dbeta(1:49*0.02, shape1=object@result$fitNULL$thetaNULL[1], shape2=object@result$fitNULL$thetaNULL[2]), col="green", lwd=4, lty=2)
			} else if(what=="fitBM") {
				if(is.null(object@result$fitBM))
					stop("PAN-error: 'fitBM' result not available yet!")
				##precompute parameters for graphics
				parts<-unique(object@partition)
				partCounts<-sapply(parts, function(part) sum(object@partition==part))
				nrPart<-length(parts)
				if(nrPart>1) {
					if(nrPart%%2==0)
						layout.series<-1:nrPart
					else 
						layout.series<-c(1:nrPart, NULL)
					layout(matrix(layout.series, length(layout.series)/2, 2))
				}
				for(s in 1:nrPart) {
					##check fitting of model to data
					hist(object@association[which(object@partition==parts[s])], freq=FALSE, breaks=50, col=gray(0.8), 
							border=gray(0.8), xlab="Transformed cosine similarity", main=paste("Beta-Mixture Model Fitting (",parts[s],")",sep=""))
					lines(1:49*0.02,sapply(1:49*0.02, betaMixture_overall_pdf, PI=object@result$fitBM$pi[s, ], 
									p=object@result$fitBM$theta),type="l", col=gray(0.3), lwd=5, lty=2)
					lines(1:49*0.02,sapply(1:49*0.02, betaMixture_partial_pdf, pi=object@result$fitBM$pi[s, 1], 
									a=object@result$fitBM$theta[1], b=object@result$fitBM$theta[2]),type="l", col="blue", lwd=4, lty=2)
					lines(1:49*0.02,sapply(1:49*0.02, betaMixture_partial_pdf, pi=object@result$fitBM$pi[s, 3], 
									a=object@result$fitBM$theta[5], b=object@result$fitBM$theta[6]),type="l", col="red", lwd=4, lty=2)
					lines(1:49*0.02,sapply(1:49*0.02, betaMixture_partial_pdf, pi=object@result$fitBM$pi[s, 2], 
									a=object@result$fitBM$theta[3], b=object@result$fitBM$theta[4]),type="l", col="green", lwd=4, lty=2)
				}
			}
		}
)
##show summary information on screen
setMethod(
		"show",
		"BetaMixture",
		function(object) {
			cat("A three-component beta-mixture model:\n")
			summarize(object, what=c("input"))
		}
)
##print summary information on screen
setMethod(
		"summarize",
		c("BetaMixture", "character_Or_missing"),
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
			if(any(c("ALL","fitNULL") %in% what)) {
				cat("-NULL fitting \n")
				if(!is.null(object@summary$fitNULL))
					print(object@summary$fitNULL, quote=FALSE, justify="centre")
				else 
					cat("not available yet!")
				cat("\n")
			}
			if(any(c("ALL","fitBM") %in% what)) {
				cat("-Beta-Mixture model fitting\n")
				if(!is.null(object@summary$fitBM)) {
					cat("--parameters: \n")
					print(object@summary$fitBM[["parameters"]], quote=FALSE, justify="centre")
					cat("--control arguments: \n")
					print(object@summary$fitBM[["control args"]], quote=FALSE, justify="centre")
					cat("--results: \n")
					print(object@summary$fitBM[["result"]], quote=FALSE, justify="centre")	
				} else 
					cat("not available yet!")
				cat("\n")
			}
		}
)







