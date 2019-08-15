######Variance optimization (VO) for optimize k selection
## Functions modified from MetaKTSP package: https://github.com/metaOmics/MetaKTSP/blob/master/R/MetaKTSP_Source.R

# expr: expression matrix (gene by sample)
# label: binary label of the expr, matched with colnames of expr
# tsp.df: dataframe output from cal.TSPscore.R
# K.max: maximum selected K


# Extract disconjuncting genes
unique.pair <- function (tsp.df, K.max=500) {
  uni.tsp <- tsp.df[1, ]
  i <- 2
  repeat {
    if (sum(tsp.df[i, c('geneX', 'geneY')] %in% c(uni.tsp$geneX, uni.tsp$geneY)) == 0) {
      uni.tsp <- rbind(uni.tsp, tsp.df[i, ])
    }
    i <- i + 1
    if (nrow(uni.tsp) == K.max) {break}
  }
  uni.tsp
}




VO.select <- function(expr, label, tsp.df, K.max, filt.pair = FALSE) {
  suppressPackageStartupMessages(library(doParallel))
  registerDoParallel()
  
  if (K.max < 2) stop('K.max should be larger than 1.')
  if (length(unique(label)) > 2) stop('Labels should only include two classes.')
  
  uni.tsp <- unique.pair(tsp.df, K.max = K.max)

  
  # Calculate t-statistics of the target function
  label <- factor(label)
  label0 <- levels(label)[1]
  label1 <- levels(label)[2]
  
  k.try <- 2:K.max
  total.tau <- foreach(k=k.try, .combine=c) %dopar% {
    Var.0 <- var(apply(foreach(ii=1:k, .combine=rbind) %do% {
        expr[uni.tsp$geneX[ii], label==label0] > expr[uni.tsp$geneY[[ii]], label==label0]
      }, 2, sum))
    
    Var.1 <- var(apply(foreach(ii=1:k, .combine=rbind) %do% {
        expr[uni.tsp$geneX[ii], label==label1] > expr[uni.tsp$geneY[[ii]], label==label1]
      }, 2, sum))
    
    sum(as.numeric(uni.tsp$score[1:k])) / sqrt(Var.0 + Var.1)
  }
  
  # optimal K
  i <- which(total.tau==max(total.tau))[1]
  print(paste0(' The optimized K = ', k.try[i]))
  if (filt.pair) {
    final.tsp <- uni.tsp[1:k.try[i], ]
  } else {
    final.tsp <- uni.tsp
  }
  
  list(opti.K = k.try[i], kpairs = final.tsp)
}

