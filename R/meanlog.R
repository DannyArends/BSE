#
# meanlog.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# De-log large numbers
#

meanlog<-function(x){
  x <- as.matrix(x)
	if(any(is.na(x))){
	  means0 <- NA
	}else{
		x1a<-as.matrix(x[x <  -700 | x  > 700])     # Why 700 ?
		x1b<-as.matrix(x[x >= -700 & x <= 700])
		if(length(x1a) == 0){
		  means0 <- log(mean(exp(x)))
		}else{
			means1 <- sum(brob(x1a))
			means2 <- length(x1b)*mean(exp(x1b))
			if(length(x1b) == 0) means0 <- log((1/length(x))*(means1))	
			if(length(x1b) >  0) means0 <- log((1/length(x))*(means1+means2))
		}
	}
	return(means0)
}
