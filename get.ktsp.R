# Function to calculate multi-class TSP scores in ways of:
# 1. One vs. all other
# 2. Pairwise

suppressPackageStartupMessages({
  library(switchBox)
})


getTSP <- function(grps, expr, labels, way) {
  tsp.list <- list()
  if (way == 'one') { # one v.s. all other
    for (grp in grps) {
      grp.vec <- rep(0, length(labels))
      grp.vec <- factor(replace(grp.vec, which(labels == grp), 1)) 
      tmp <- SWAP.CalculateScores(as.matrix(expr), grp.vec, FilterFunc = NULL)
      tmp <- data.frame(score = sort(tmp$score, decreasing = T))
      tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
      tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
      rownames(tmp) <- NULL
      tsp.list[[as.character(grp)]] <- tmp
    }
  } else if (way == 'pairwise') { # pairwise
    for (i in 1:(length(grps)-1)) {
      for (j in (i+1):length(grps)) {
        exp1 <- expr[ ,labels == as.character(grps[i])]
        exp2 <- expr[ ,labels == as.character(grps[j])]
        grp.data <- cbind.data.frame(exp1, exp2)
        grp.vec <- factor(labels[labels == as.character(grps[i]) | labels == as.character(grps[j])])
        tmp <- SWAP.CalculateScores(as.matrix(grp.data), grp.vec, FilterFunc = NULL)
        tmp <- data.frame(score = sort(tmp$score, decreasing = T))
        tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
        tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
        rownames(tmp) <- NULL
        tsp.list[[sprintf('%s_%s', grps[i], grps[j])]] <- tmp
      }
    }
  } else {print('The options of way are just "one" or "pairwise"!')}
  tsp.list
}



# Merge k gene pairs from each group
mergePairs <- function(Ks, tsp.list, no.model) {
  # Function to eliminate duplicated rows
  delete.dup <- function(table) {
    dup.table <- table[which(duplicated(table[, c(1, 2)])), ]
    if (nrow(dup.table) != 0) {
      sum.dup <- unique(dup.table %>% group_by_(.dots=c("geneX", "geneY")) %>% 
                          mutate(Dclass=paste(class, collapse = ",")) %>% select(-class))
      table <- table[which(!duplicated(table[, c(1, 2)])), ]
      table$class <- as.vector(table$class)
      for (n in 1:nrow(sum.dup)) {
        row <- sum.dup[n, ]
        rep.names <- strsplit(row$Dclass, ",")[[1]]
        ori.name <- as.character(table[(table$geneX==row$geneX & table$geneY==row$geneY),]$class)
        if (!ori.name %in% rep.names) {
          ori.name <- paste(c(ori.name, rep.names), collapse = ",")
        } else {
          ori.name <- paste(rep.names, collapse = ",")
        }
        table[(table$geneX==row$geneX & table$geneY==row$geneY),]$class <- ori.name
      }
    }
    table
  }
  
  kpairs <- list()
  # Extract k pairs
  for (K in Ks) {
    ifelse (no.model == 3, kpairs[[as.character(K)]] <- list(), kpairs[[as.character(K)]] <- data.frame())
    for (s in 1:length(tsp.list)) {
      name <- names(tsp.list)[s]
      scores <- tsp.list[[s]]
      if (no.model == 3) { # For model 3
        kpairs[[as.character(K)]][[name]] <- scores[1:K, ] %>% mutate(class = rep(name, K)) %>% select(-score)
      } else { # For model 1&2 merging
        table <- scores[1:K, ] %>% mutate(class = rep(name, K)) %>% select(-score)
        kpairs[[as.character(K)]] <- rbind.data.frame(kpairs[[as.character(K)]], table)
      }
    }
  }
  
  # Delete repeat rows
  if (no.model != 3) {
    for (i in 1:length(kpairs)) {
      kpairs[[i]] <- delete.dup(kpairs[[i]])
    }
  }
  kpairs
}
