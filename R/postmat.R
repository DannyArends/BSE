#
# postmat.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# Read Metropolis output
#

postmat<-function(eachtemp, s, k){
  draws<-eachtemp$sims.matrix 
	Rsample<-eachtemp$n.keep

	amatrix<-matrix(draws[,1:s], Rsample, s)

	bmatrix<-array(NA, dim=c(Rsample, s, k))
	bmatrix[, , 1]<-draws[, (s+1):(2*s)]
	bmatrix[, 1, 2]<-0
	bmatrix[, 2:s, 2]<-draws[, (2*s+1):(3*s-1)]

	return(list(amatrix=amatrix, bmatrix=bmatrix))
}
