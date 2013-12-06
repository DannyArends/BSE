#
# BSE.example.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# Run example
#

BSE.example <- function(){
  y <-as.matrix(read.table("input_66002.csv"))
	s<-ncol(y)
	g<-nrow(y)
	k<-2   #number of hypothesed latent variables  
	npoints<-30

	R     <-  600 	# n.iter
 	Rburn <-  100   # n.burnin 
	Rthin <-  1		  # n.thin  

	# Metropolis 
	Metro.draws<-Metropolis(y,k,R,Rburn,Rthin) 
	post.out<-postmat(Metro.draws,s,k)
	amatrix<-post.out$amatrix
	bmatrix<-post.out$bmatrix
	
	estimators<-margGH.bs(y,k,amatrix,bmatrix,npoints) 
}