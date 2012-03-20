##cosine similarity
cosineSim <- function(x) {
	nom <- t(x) %*% x
	den <- sqrt(apply(x^2,2,sum))
	cos <- nom / den %*% t(den) 
	##return(as.dist(1-cos))
	return(cos)
}
##cosine distance
cosineDist<-function(x) {
	as.dist(1-cosineSim(x))
}
##compute association scores
###... would be used for the arguments of function 'cor' 
assoScore<-function(pheno, metric="cosine", upperTri=TRUE, transform=TRUE, verbose=TRUE, ...) {
			assoScoreSingle<-function(pheno, ...) {
				##compute association scores
				if(metric=="cosine") {
					asso<-cosineSim(pheno)
				} else if(metric=="correlation") {
					asso<-cor(pheno, ...)
				} else {
					stop("PAN-error: only cosine and correlation coefficients supported by current version")
				}
				asso
			}

			asso<-assoScoreSingle(pheno, ...)
			##avoid R exceeding association ranges [-1, 1]
			asso[asso<(-1)]<-(-1)
			asso[asso>1]<-1
			##transform to [0, 1]
			if(transform)
				asso<-(asso+1)/2
			##return
			if(upperTri)
				return(asso[upper.tri(asso)])
			else
				return(asso)
		}


