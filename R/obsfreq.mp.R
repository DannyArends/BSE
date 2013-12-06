#
# obsfreq.mp.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# Get observed patterns
#

obsfreq.mp <- function(x){
  s<-ncol(x)
	pats <- apply(x, 1, paste, collapse = "")
	freqs <- table(pats)
	obs.mp <-as.matrix(freqs)
	n.mp<-dim(obs.mp)[1]
	i.obsm <- apply(cbind(dimnames(obs.mp)[[1]]), 1, function(x) {
                    nx<-nchar(x)
 					          out <- substring(x, 1:nx, 1:nx)
   					        out
   					 }) 
	pm.mp<- matrix(as.numeric(t(i.obsm)),ncol=s)
	return(list(obs.mp=obs.mp, n.mp=n.mp, pm.mp=pm.mp))
}
