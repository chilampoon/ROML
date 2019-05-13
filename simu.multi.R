# The function to generate simulation data

suppressPackageStartupMessages({
  library(truncnorm)
})


simu.multi <- function(g.null, g.sig, alpha, n.train, n.test, SDnullMu, DFnullSD, lb.shift, ub.shift, testMu.shift, testSD.shift, pShift, seed=15213) {
  ### input
  # g.null - number of null (non-various) genes
  # g.sig - number of significant (various) genes
  # alpha - rate of only-one-class-DE-gene (1, 0, 0); while 1-alpha is for the rate of all-diff-DE-gene (-1, 0, 1)
  # n.train - a vector for number of training samples for each class (can be multi-class)
  # n.test - a vector for number of testing samples (same length with n.train)
  # SDnullMu - SD for null gene mu value
  # DFnullSD - degree of free for null gene SD
  # lb.shift - significant gene shift lower bound
  # ub.shift - significant gene shift upper bound
  # testMu.shift - test shift mean value
  # testSD.shift - test shift SD value
  # seed - seed for random number
  
  ### output
  # Xtrain - training expression table, gene by sample
  # Ytrain - training outcome label
  # Xtest - testing expression table, gene by sample
  # Ytest - testing outcome label
  # sig.idx - significant gene index, TRUE for sig and FALSE for non-sig
  # geneMUtrain - gene mu value, a gene-by-group matrix, for training data
  # geneMUtest - gene mu value, a gene-by-group matrix, for testing data because of the shift
  # geneSD - gene SD value, a vector with length of gene
  # geneShift - TRUE or FALSE showing whether this is a shift gene
  
  set.seed(seed)
  
  g.sig <- as.integer(g.sig)
  if (g.sig <= 0) stop("g.sig - number of significant genes - has to be a positive integer")
  
  g.null <- as.integer(g.null)
  if (g.null <= 0) stop("g.null - number of null genes - has to be a positive integer")
  
  if (alpha < 0 | alpha > 1) stop("Ratio alpha has to be [0, 1]")
  
  n.train <- as.integer(n.train)
  if (length(n.train) != 3) stop("Training data have to have 3 outcomes, that is length(Ntrain) = 3")
  
  n.test <- as.integer(n.test)
  if (length(n.train) != length(n.test)) stop("Length of n.train has to be equal to length of n.test.")
  
  
  # prepare gene number and index
  g.all <- g.sig + g.null
  g.total <- sum(g.all) #??
  g.start.idx <- c(1, g.sig+1)
  g.end.idx <- c(g.sig, g.all)
  
  
  # prepare sample number and index
  n.class <- length(n.train)  # number of outcome categories
  n.all <- n.train + n.test # total number
  n.total <- sum(n.all) 
  n.start.idx <- rep(1, times=length(n.all)) # start index
  n.end.idx <- rep(n.all[1], times=length(n.all)) # end index
  for (i in 2:length(n.all)) {
    n.start.idx[i] <- n.start.idx[i-1] + n.all[i-1] ## ??
    n.end.idx[i] <- n.end.idx[i-1] + n.all[i] ## ??
  }
  
  
  ## 1. Simulate null/baseline genes (exprssion not vary among all classes)
  E <- matrix(0, nrow=g.total, ncol=n.total)  # expression table
  gene_mu <- rnorm(g.total, mean=0, sd=SDnullMu)  # gene mu follows normal distribution
  geneSD <- sqrt(rchisq(n=g.total, df=DFnullSD))   # gene variance follows chisq distribution
  for (i in 1:g.total) {
    E[i, ] <- rnorm(n.total, mean=gene_mu[i], sd=geneSD[i]) 
  }
  
  
  ## 2. Simulate significant genes (vary in 1 - one vs all other or 2 - pairwise ways) 
  # 2.1 only-one-DE gene, alpha type, type A
  g.alpha <- as.integer(g.sig * alpha) # number of one-other-sig-genes
  if (g.alpha > 0) {
    shiftC1 <- matrix(rtruncnorm(g.alpha*n.class, a=lb.shift, b=ub.shift, mean=0, sd=SDnullMu), nrow=g.alpha, ncol=n.class)
    
    directionCandidate1 <- rbind(c(1,0,0), c(0,1,0), c(0,0,1), c(-1,0,0), c(0,-1,0), c(0,0,-1))
    direction1 <- matrix(0, nrow=g.alpha, ncol=n.class)  # showing the gene direction for each outcome category
    for (i in 1:g.alpha) {
      direction1[i, ] <- directionCandidate1[sample(nrow(directionCandidate1), 1), ] # randomly choose from candidate
      
      for (j in 1:n.class) {
        E[i, n.start.idx[j]:n.end.idx[j]] <- E[i, n.start.idx[j]:n.end.idx[j]] + direction1[i, j] * shiftC1[i, j]
      }
    }
  } else {
    direction1 <- c()
    shiftC1 <- c()
  }
  
  # 2.2 all-DE gene, 1-alpha type, type B
  g.alphaN <- g.sig - g.alpha # number of pairwise-sig-genes
  if (g.alphaN > 0) {
    shiftC2 <- matrix(rtruncnorm(g.alphaN*n.class, a=lb.shift, b=ub.shift, mean=0, sd=SDnullMu), nrow=g.alphaN, ncol=n.class)
    directionCandidate2 <- rbind(c(1,0,-1), c(1,-1,0), c(0,-1,1), c(0,1,-1), c(-1,1,0), c(-1,0,1))
    direction2 <- matrix(0, nrow=g.alphaN, ncol=n.class)  # showing the gene direction for each outcome category
    for (i in 1:g.alphaN) {
      direction2[i, ] <- directionCandidate2[sample(nrow(directionCandidate2), 1), ] # randomly choose from candidate
      
      for (j in 1:n.class) {
        E[g.alpha+i, n.start.idx[j]:n.end.idx[j]] <- E[i+g.alpha, n.start.idx[j]:n.end.idx[j]] + direction2[i,j] * shiftC2[i,j]
      }
    }
  } else {
    direction2 <- c()
    shiftC2 <- c()
  }
  
  
  ## summarize the data
  # sig gene index
  sig.idx <- rep('N', times=g.total)
  sig.idx[1:g.sig] <- rep(c('A', 'B'), times=c(g.alpha, g.alphaN))
  
  # gene effect
  direction <- rbind(direction1, direction2)
  shiftC <- rbind(shiftC1, shiftC2)
  geneMUtrain <- matrix(gene_mu, nrow=g.total, ncol=n.class)
  if (g.sig > 0) {
    geneMUtrain[1:g.sig, ] <- geneMUtrain[1:g.sig, ] + direction * shiftC
  }  
  geneMUtest <- geneMUtrain
  
  ## split training and testing expression
  train.idx <- c()
  for (j in 1:n.class) {
    train.idx <- c(train.idx, (n.start.idx[j]:(n.start.idx[j ] + n.train[j] - 1)))
  }
  
  ifelse(length(train.idx) > 0, Xtrain <- E[ ,train.idx], Xtrain <- c())
  
  # Select shift index, p% genes in test data will add shift.
  shift.idx <- sample(1:g.total, round(g.total * pShift))
  gene.shift <- rep(FALSE, times=g.total)
  gene.shift[shift.idx] <- TRUE
  
  # add shift to test value
  if (length(train.idx) < n.total) {
    Xtest <- E[, -train.idx]
    for (i in shift.idx) {
      shiftTmp <- rnorm(1, mean=testMu.shift, sd=testSD.shift)
      Xtest[i, ] <- Xtest[i, ] + shiftTmp
      geneMUtest[i, ] <- geneMUtest[i,] + shiftTmp
    }
  } else {
    Xtest <- c()
  }
  
  Xtrain <- as.data.frame(Xtrain)
  Xtest <- as.data.frame(Xtest)
  ###### the rowname set alpha and beta
  rownames(Xtrain) <- sapply(seq_along(rownames(Xtrain)), function(x) paste0(sig.idx[x], x)) -> rownames(Xtest)
  rownames(geneMUtrain) <- rownames(Xtrain)
  colnames(Xtrain) <- c(1:ncol(Xtrain))
  colnames(Xtest) <- c(1:ncol(Xtest))
  
  # train and test y label
  Ytrain <- rep(0:(n.class-1), times=n.train)
  Ytest <- rep(0:(n.class-1), times=n.test)
  
  list(Xtrain=Xtrain, Ytrain=factor(Ytrain), Xtest=Xtest, Ytest=factor(Ytest), geneMUtrain=geneMUtrain, geneMUtest=geneMUtest,
       sig.idx=sig.idx, geneSD=geneSD, geneShift=gene.shift)
}