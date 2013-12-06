#
# gh.each.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# Approximate each integral using G-H
#

gh.each<-function(x,a,b,k,npoints){
  k<-k
	s <- dim(x)[2]
	g <- dim(x)[1]
	obspat<-obsfreq.mp(x)
	pm<-obspat$pm.mp
	n<-obspat$n.mp
	times<-obspat$obs.mp

	a <- matrix(data = a, nrow = s, ncol = 1)
	b <- matrix(data = b, nrow = k, ncol = s)
	n <- n

	out <- gauss.quad.prob(npoints, "normal") #### NZ = Number of points
	zpoints<-out$nodes
	zweights<-out$weights

	zv <- matrix(data = zpoints, nrow = npoints, ncol = 1)
	wz <- matrix(data = zweights, nrow = npoints, ncol = 1)
	np<-length(zv)
			
	bar <- array(b,dim=c(1,k,s))
	bar2 <- array(NA,dim=c(np^k,k,s))
	for(i in 1:s){
	  bar2[,,i]<-rep(bar[,,i], each = np^k)
	}

	zvma2 <- array(zv,dim=c(np,1,1))
	wzmat <- matrix(NA,np^k,k)
			

	zmat <- matrix(NA,np^k,k)
	rob  <- matrix(nrow = n, ncol = 1)

	for (i in 1:k){
    zmat[,1]<-rep(zv ,np^(k-1))
		zmat[,i]<-rep(zv ,each=np^(i-1))

		wzmat[,1]<-rep(wz,np^(k-1))
		wzmat[,i]<-rep(wz,each=np^(i-1))
  }
  zarray <- array(zmat,c(np^k,k,s))
	byit   <- bar2 * zarray
	byit2  <- rowSums(aperm(byit,c(1,3,2)), na.rm = FALSE, dims = 2)
	wzm<-as.matrix(apply(wzmat,1,prod))

	azb<-t(matrix(a, s,np^k )) + byit2
	e <- exp(azb)
	pk <-  e / (1 + e)

  # TODO: Issue a warning ?
	pk[pk>=0.999999] = 0.9999999
	pk[pk<=0] = 0.00000000000000000001
	
	ptrn2 <- array(rep(t(pm),each=np^k),dim=c(np^k,s,n))
	ps23  <- array(rep(pk),dim=c(np^k,s,n))

	gh23   <- colSums(aperm(ptrn2*log(ps23)+(1-ptrn2)*log(1-ps23),c(2,1,3)), na.rm = FALSE, dims = 1)
	abzw23 <- gh23 + log(wzm[ , rep(seq_len(ncol(wzm)), n)])

 	prob <- as.matrix(apply(exp(abzw23), 2, sum))

	return(sum(log(prob)*times))
}
