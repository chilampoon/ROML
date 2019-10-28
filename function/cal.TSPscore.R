# Calculate TSP scores 
# expr: gene-by-sample
source('~/ktsp/script/function/pair.t.stat.R')

## add multiclass model calculation method

cal.TSPscore <- function(expr, label, t.stat=TRUE, multi.model=FALSE, K.max=10000, save=FALSE, path=NULL) {
  suppressPackageStartupMessages(library(switchBox))
  grp <- factor(label)
  if (multi.model == FALSE) { 
    # 1. binary & model3
    tsp <- SWAP.CalculateScores(as.matrix(expr), grp, FilterFunc = NULL)
    tsp <- data.frame(score = sort(tmp$score, decreasing = T))
    tsp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
    tsp$geneY <- gsub('^(.*?),', '', rownames(tmp))
    rownames(tsp) <- NULL
    tsp.list <- tsp[tsp$score > 0, ] # use positive scores only, since they are symmetric
  } else if (multi.model == '1') {
    # 2. one-vs-other
    tsp.list <- list()
    for (g in grp) {
      grp.vec <- rep(0, length(label))
      grp.vec <- factor(replace(grp.vec, which(label == g), 1)) 
      tsp <- cal.TSPscore(expr, grp.vec)
      tsp.list[[g]] <- tsp[1:K.max,]
      if (save == TRUE) save(tsp.list, file = path)
    }
  } else if (multi.model == '2') {
    # 3. pairwise
    # split expr data
    expr <- cbind.data.frame(t(expr), type=grp)
    sub.expr <- list()
    for (g in grp) {
      sub.expr[[g]] <- expr[expr$type == g, ]
    }
    
    tsp.list <- list()
    for (i in 1:length(grp)) {
      for (j in (i+1):(length(grp)-1)) {
        grp.data <- rbind.data.frame(sub.expr[[grp[i]]], sub.brca.expr[[grp[j]]])
        grp.vec <- factor(grp.data$type)
        grp.data <- t(as.matrix(grp.data[, which(colnames(grp.data) != "type")]))
        tsp <- cal.TSPscore(grp.data, grp.vec)
        tsp.list[[paste0(grp.names[i], "_", grp.names[j])]] <- tsp[1:K.max, ]
        if (save == TRUE) save(tsp.list, file = path)
      }
    }
  }
  
  
  # if (t.stat == TRUE) {
  #   t.df <- as.data.frame(t(apply(tsp[ ,-1], 1, split_statistic, data.train=expr, label=label)))
  #   tsp <- cbind(tsp, t.df)
  # }
  tsp.list
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

