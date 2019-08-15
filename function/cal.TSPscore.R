# Calculate TSP scores 

# Here only for binary outcome

cal.TSPscore <- function(expr, labels) {
  suppressPackageStartupMessages(library(switchBox))
  grp.vec <- factor(labels)
  tmp <- SWAP.CalculateScores(as.matrix(expr), grp.vec, FilterFunc = NULL)
  tmp <- data.frame(score = sort(tmp$score, decreasing = T))
  tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
  tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
  rownames(tmp) <- NULL
  tmp <- tmp[tmp$score > 0, ] # use positive scores only, since they are symmetric
  tmp
}

