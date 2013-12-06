#
# BSE.example.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# k: Number of hypothesed latent variables  
#

BSE.example <- function(k = 2, n.points = 30, n.iter = 600, n.burnin = 100, n.thin = 1){
  data(input66002)
  s <- ncol(input66002)

  # Call the metropolis function wrapper
  metro.draws <- metropolis(input66002, k, n.iter, n.burnin, n.thin) 

  post.out <- postmat(metro.draws, ncol(input66002), k)

  estimators<-margGH.bs(y, k, post.out$amatrix, post.out$bmatrix, n.points) 

  cat("Predicted: ", estimators, "\n", sep="")
}