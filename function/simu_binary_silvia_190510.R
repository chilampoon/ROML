### data simulation for TSP 
### simulate the data follow a certain distribution, so that truth is known

library(truncnorm)

simu.multi <- function(Gnull, Gsig, Ntrain, Ntest, SDgenemu=5, DFgenesd=1, Clow=2, Cup=Inf, MuShift=0, SDshift=1, pShift=0.3,seed=15213){
  ### input
  # Gnull - number of null genes
  # Gsig - number of significant genes
  # Ntrain - a vector for number of training samples for each class (can be multi-class)
  # Ntest - a vector for number of testing samples (same length with Ntrain)
  # SDgenemu - SD for null gene mu value
  # DFgenesd - degree of free for null gene SD
  # Clow - significant gene shift lower bound
  # Cup - significant gene shift upper bound
  # MuShift - test shift mean value
  # SDshift - test shift SD value
  # pShift - proportion of the genes to shift in test
  # seed - seed for random number
  
  ### output
  # Xtrain - training expression table, gene by sample
  # Ytrain - training outcome label
  # Xtest - testing expression table, gene by sample
  # Ytest - testing outcome label
  # sigIndex - significant gene index, TRUE for sig and FALSE for non-sig
  # geneMUtrain - gene mu value, a gene-by-group matrix, for training data
  # geneMUtest - gene mu value, a gene-by-group matrix, for testing data because of the shift
  # geneSD - gene SD value, a vector with length of gene
  # geneShift - TRUE or FALSE showing whether this is a shift gene
  
  set.seed(seed)
  
  Gsig=as.integer(Gsig)
  if(Gsig<=0){
    stop("Gsig - number of significant genes - has to be a positive integer")
  }
  
  Gnull=as.integer(Gnull)
  if(Gnull<=0){
    stop("Gnull - number of null genes - has to be a positive integer")
  }
  
  Ntrain=as.integer(Ntrain)
  if(length(Ntrain)<2){
    stop("Training samples need to have at least two outcomes.")
  }
  
  Ntest=as.integer(Ntest)
  if(length(Ntrain)!=length(Ntest)){
    stop("Length of Ntrain has to be equal to length of Ntest.")
  }
  
  if(pShift<0 | pShift>1){
    stop("pShift is a proportion in [0,1]")
  }
  
  ### prepare gene number and index
  Gall=Gsig+Gnull
  Gtotal=sum(Gall)
  GstartInd=c(1,Gsig+1)
  GendInd=c(Gsig,Gall)
  
  ### prepare sample number and index
  K=length(Ntrain)  # number of outcome categories
  Nall=Ntrain+Ntest # total number
  Ntotal=sum(Nall)
  NstartInd=rep(1,times=length(Nall)) # start index  
  NendInd=rep(Nall[1],times=length(Nall)) # end index
  for(i in 2:length(Nall)){
    NstartInd[i]=NendInd[i-1]+1
    NendInd[i]=NstartInd[i]+Nall[i]-1
  }
  
  ### null distribution
  E=matrix(0,nrow=Gtotal,ncol=Ntotal)  # expression table
  gene_mu=rnorm(Gtotal,mean=0,sd=SDgenemu)  # gene mu follows normal distribution
  gene_sd=sqrt(rchisq(n=Gtotal,df=DFgenesd))   # gene variance follows chisq distribution
  for(i in 1:Gtotal){
    E[i,]=rnorm(Ntotal,mean=gene_mu[i],sd=gene_sd[i]) 
  }
  
  ### sig gene
  shiftC=matrix(rtruncnorm(Gsig*K,a=Clow,b=Cup,mean=0,sd=SDgenemu),nrow=Gsig,ncol=K)
  direction=matrix(0,nrow=Gsig,ncol=K)  # showing the gene direction for each outcome category
  for(i in 1:Gsig){
    while(length(unique(direction[i,]))==1){  # not the sample direction
      direction[i,]=sample(c(-1,0,1),size=K,replace=T)
    }
    
    for(j in 1:K){
      E[i,NstartInd[j]:NendInd[j]]=E[i,NstartInd[j]:NendInd[j]] + direction[i,j]*shiftC[i,j]
    }
  }
  
  ## summarize the data
  
  # sig gene index
  sigIndex=rep(FALSE,times=Gtotal)
  sigIndex[1:Gsig]=TRUE
  
  # gene effect
  geneMUtrain=matrix(gene_mu,nrow=Gtotal,ncol=K)
  geneMUtrain[sigIndex,]=geneMUtrain[sigIndex,]+direction*shiftC
  geneMUtest=geneMUtrain
  
  ## split training and testing expression
  trainInd=c()
  for(j in 1:K){
    trainInd=c(trainInd,(NstartInd[j]:(NstartInd[j]+Ntrain[j]-1)))
  }
  
  if(length(trainInd)>0){
    Xtrain=E[,trainInd]
  }else{
    Xtrain=c()
  }
  
  # select shift Index
  shiftInd=sample(1:Gtotal,round(Gtotal*pShift))
  geneShift=rep(FALSE,times=Gtotal)
  geneShift[shiftInd]=TRUE
  
  if(length(trainInd)<Ntotal){
    Xtest=E[,-trainInd]
    ## add shift to test value
    for(i in shiftInd){
      shiftTmp=rnorm(1,mean=MuShift,sd=SDshift)
      Xtest[i,]=Xtest[i,]+shiftTmp
      geneMUtest[i,]=geneMUtest[i,]+shiftTmp
    }
  }else{
    Xtest=c()
  }
  
  ## train and test y label
  Ytrain=rep(1:K,times=Ntrain)
  Ytest=rep(1:K,times=Ntest)
  
  
  return(list(Xtrain=Xtrain,Ytrain=Ytrain,Xtest=Xtest,Ytest=Ytest,
              sigIndex=sigIndex, geneMUtrain=geneMUtrain, geneMUtest=geneMUtest,
              geneSD=gene_sd, geneShift=geneShift))
  
}



#### testing the simu.multi function, very small number
Gnull=50
Gsig=10
Ntrain=c(10,10)
Ntest=c(5,10)
# SDgenemu=5
# DFgenesd=1
# Clow=2
# Cup=Inf
# MuShift=0
# SDshift=1
# seed=15213

out=simu.multi(Gnull, Gsig, Ntrain, Ntest,pShift=0.5, SDshift=5)

# ### plot to check gene expression
E=out$Xtest
E=out$Xtrain
plot(NA,xlim=c(1,ncol(E)),ylim=range(E),xlab="sample",ylab="expression")
for(i in (Gsig+1):nrow(E)){
  lines(x=1:ncol(E),y=E[i,],col="grey")
}
for(i in 1:Gsig){
    lines(x=1:ncol(E),y=E[i,],col=sample(rainbow(10),1))
}
