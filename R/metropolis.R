#
# metropolis.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# Parameter estimation: fit the model in winbugs
#

metropolis<-function(y, k, n.iter = 600, n.burnin = 100, n.thin = 1){
  g<-nrow(y)		# Number of individuals
  s<-ncol(y)		# Number of items (observed variables)	
	
	# Initial values
	ini_a  <- c(-1.380, -1.073, -0.549,  0.439,  1.008,  2.425)
	ini_b1 <- c( 1.821,  1.944,  1.346, -1.278, -1.581, -2.342)
	ini_b2 <- c(    NA,  1.034,  1.167,  3.072,  1.538,  2.067)  # NA as value seems weird

	ini_z1 <- rnorm(g)
	ini_z2 <- rnorm(g)
		 
	parameters <- c("a","b1","b2")	  # Parameters to the model ?
	datas      <- list("g","s","y")   # Temperament is here, variables in the model ?
	inits <- function(){
	  list(a=ini_a, b1=ini_b1, b2=ini_b2, z1=ini_z1, z2=ini_z1)
	} #changes at each T

	eachtemp <- bugs(datas, inits, parameters, "2factor.bugs", n.chains = 1, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin)
	return(eachtemp)
}
