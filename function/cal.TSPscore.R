# Calculate TSP scores 
# Here only for binary outcome
source('~/ktsp/script/function/pair.t.stat.R')


cal.TSPscore <- function(expr, label, t.stat=TRUE) {
  suppressPackageStartupMessages(library(switchBox))
  grp.vec <- factor(label)
  tmp <- SWAP.CalculateScores(as.matrix(expr), grp.vec, FilterFunc = NULL)
  tmp <- data.frame(score = sort(tmp$score, decreasing = T))
  tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
  tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
  rownames(tmp) <- NULL
  tmp <- tmp[tmp$score > 0, ] # use positive scores only, since they are symmetric
  
  if (t.stat == TRUE) {
    t.df <- as.data.frame(t(apply(tmp[ ,-1], 1, split_statistic, data.train=expr, label=label)))
    tmp <- cbind(tmp, t.df)
  }
  tmp
}


# Extract disconjuncting genes
unique.pair <- function (tsp.df, K.max=5000) {
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

