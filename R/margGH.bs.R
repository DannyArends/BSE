#
# margGH.bs.R
#
# copyright (c) 2013-2015 - Silia Vitoratou and Danny Arends
# last modified Dec, 2013
# first written Dec, 2013
# 
# The function margGH.bs takes the posterior samples from the function metropolis
# and a) creates the reference function g(theta) and b) computes the BS marginal 
# likelihood estimators.
#

margGH.bs <- function(y,k,amatrix,bmatrix,npoints){
  R<-nrow(amatrix)
  L<-R
  g <- nrow(y)		# Number of individuals
  s <- ncol(y)		# Number of items (observed variables)

  #-----------------------------------------------------------------------------#
  #        L-draws: sampling from the envelope function   g(theta)              #
  #-----------------------------------------------------------------------------#
  # Draw parameters from posterior sample

  meana <- apply(amatrix,2,mean)
  aLmatrix <- rmvnorm(L,meana,cov(amatrix))

  # Break bmatrix into 2 pieces: normal and lognormal values

  btrunk<-matrix(NA,dim(bmatrix)[1],dim(bmatrix)[3])
  for (f in 1:k){ btrunk[,f]<-bmatrix[,f,f] }

  lbtrunk <- log(btrunk)
  bnormal <- bmatrix[,2:s,1]

  if(k>1){
    for(f in 2:k){
      m<-f+1
      step<-bmatrix[,m:s,f]
      bnormal<-cbind(bnormal,step)
    }
  }

  bMCatrixT<-btrunk
  bMCatrixN<-bnormal


  # Generate bLmatrixm from g Normal

  bmatrixm  <- cbind(bnormal,lbtrunk)
  meanb     <- apply(bmatrixm,2,mean)
  covbpost  <- cov(bmatrixm)
  bLmatrixm <- rmvnorm(L,meanb,covbpost)


  # Have bLmatrix into 2 peaces for prior and proposal
  if(k==1){
    lbLtrunk<-matrix(bLmatrixm[,(ncol(bLmatrixm)-k+1):ncol(bLmatrixm)])
  }else if(k>1){ 
    lbLtrunk<-bLmatrixm[,(ncol(bLmatrixm)-k+1):ncol(bLmatrixm)]
  }
  bLnormal  <- bLmatrixm[,1:(ncol(bLmatrixm)-k)]
  bLmatrixT <- lbLtrunk
  bLmatrixN <- bLnormal

  # Structure bLmatrix for likelihood (exp bLtrunk first)

  bLtrunk<-exp(lbLtrunk)		# bLtrunk follows lognormal

  bLmatrix<-array(NA,dim=c(L,s,k))
  u <- s-1
  l <- 1
  for (f in 1:k){
    bLmatrix[,f,f]<-bLtrunk[,f]
    bLmatrix[,(f+1):s,f]<-bLnormal[,l:u]

    l <- u+1
    u <- u+s-f-1
  }

  if(k>1){ 
    for(i in 2:k){
      bLmatrix[,1:(i-1),i] <- 0
    }
  }

  ### notes ###
  # bmatrixm has logbetas-hence all normal
  # bLmatrixm has logbetas-hence all normal
  # bLmatrix, bMCatrix have normal and lognormal parameters
  # btrunk follows lognormal-lbtrunk follows normal
  # bLtrunk follows lognormal-lbLtrunk follows normal


  #-----------------------------------------------------------------------------#
  #       Compute components of BS estimators				                            #
  #-----------------------------------------------------------------------------#

  #-----------------------------------------------------------------------------#
  # G-density L draws: calculate the density on log scale                       #
  #-----------------------------------------------------------------------------#

  # Zero column not included
  lga.L <- dmvnorm(aLmatrix,meana,cov(amatrix),log=TRUE)
  lgb.L <- as.matrix(dmvnorm(bLmatrixm,meanb,covbpost,log=TRUE))

  lgtheta.L <- lga.L+lgb.L
  write(t(lgtheta.L), file = paste("",s,"_",g,"_",k,"_R =",R,"_lgtheta.L.R"), ncolumns = 1, append = TRUE, sep = " ")

  #------------------------------------------------------------------------------#
  #Prior density of L draws: put the L draws on typical priors and find their density#


  lpa.L  <- as.matrix(dmvnorm(aLmatrix,rep(0, s),4*diag(s),log=TRUE))
  lpb.LN <- as.matrix(dmvnorm(bLmatrixN,rep(0,  ncol(bLmatrixN)),4*diag( ncol(bLmatrixN)),log=TRUE))
  lpb.LT <- as.matrix(apply(log(dlnorm(bLtrunk, 0,1)),1,sum))

  lprior.L <- lpa.L+lpb.LN+lpb.LT


  write(t(lprior.L), file = paste("",s,"_",g,"_",k,"_R =",R,"_lprior.L.R"),ncolumns = 1, append = TRUE, sep = " ")

  #G-density of MCMC draws#

  lga.MC <- as.matrix(dmvnorm(amatrix,meana,cov(amatrix),log=TRUE))
  lgb.MC <- as.matrix(dmvnorm(bmatrixm,meanb,covbpost,log=TRUE))

  lgtheta.MC <- lga.MC+lgb.MC

  write(t(lgtheta.MC), file = paste("",s,"_",g,"_",k,"_R =",R,"_lgtheta.MC.R"),ncolumns = 1, append = TRUE, sep = " ")

  #------------------------------------------------------------------------------#
  #Prior density of MCMC draws#

  lpa.MC<-as.matrix(dmvnorm(amatrix,rep(0, s),4*diag(s),log=TRUE))
  lpb.MCN<-as.matrix(dmvnorm(bMCatrixN,rep(0, ncol(bMCatrixN)),4*diag( ncol(bMCatrixN)),log=TRUE))
  lpb.MCT<-as.matrix(apply(log(dlnorm(btrunk, 0,1)),1,sum))

  lprior.MC <- lpa.MC+lpb.MCN+lpb.MCT
  timegh2<-proc.time()

  write(t(lprior.MC), file = paste("",s,"_",g,"_",k,"_R =",R,"_lprior.MC.R"),ncolumns = 1, append = TRUE, sep = " ")


  #Likelihood of L and MCMC draws#
  loglike.L <- gh.all(y,aLmatrix,bLmatrix,k,npoints)
  write(t(loglike.L), file = paste("",s,"_",g,"_",k,"_R =",R,"_loglike.L.R"),ncolumns = 1, append = TRUE, sep = " ")	

  loglike.MC <- gh.all(y,amatrix,bmatrix,k,npoints)
  write(t(loglike.MC), file = paste("",s,"_",g,"_",k,"_R =",R,"_loglike.MC.R"),ncolumns = 1, append = TRUE, sep = " ")

  #-----------------------------------------------------------#
  #       Compute BS estimators                               #
  #-----------------------------------------------------------#

  HM     <- (-meanlog(as.matrix(-loglike.MC)))
  RM     <- (-meanlog(lgtheta.MC-loglike.MC-lprior.MC))
  numg   <- meanlog(0.5*(loglike.L+lprior.L-lgtheta.L))
  denomg <- meanlog(-0.5*(loglike.MC+lprior.MC-lgtheta.MC))
  BG     <- numg-denomg

  numhar   <- meanlog(-lgtheta.L)
  denomhar <- meanlog(-loglike.MC-lprior.MC)
  BH<-(numhar-denomhar)	

  results<-matrix(NA,4,1)
  results[1,1]<-HM			 	
  results[2,1]<-RM			 
  results[3,1]<-BH			 
  results[4,1]<-BG

  colnames(results)<-c("Estimated BML")
  rownames(results)<-c("HM","RM","BH","BG")
  results<-round(results,3)

  return("Bridge family estimators"=results)
}
