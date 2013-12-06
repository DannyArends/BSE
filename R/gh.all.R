#
# gh.all.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# Get likelihood for total number of observed patterns
#

gh.all<-function(x,amatrix,bmatrix,k,npoints){
  R   <- dim(bmatrix)[1]
  s   <- dim(bmatrix)[2]
  ret <- matrix(NA, R, 1)

  for(r in 1:R){
    a <- amatrix[r,]
    b <- t(bmatrix[r,,])
    ret[r,1] <- gh.each(x, a, b, k, npoints)
  }

  return(ret)
}
