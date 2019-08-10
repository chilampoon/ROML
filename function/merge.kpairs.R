# Subset k pairs from top score pairsm, input a vector of k

# Function to eliminate duplicated rows
delete.dup <- function(table) {
  dup.table <- table[which(duplicated(table[ ,c(1, 2)])), ]
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




# Merge k gene pairs from each group
merge.kpairs <- function(Ks, tsp.list, no.model) {
  kpairs <- list()

  for (K in Ks) {
    if (missing(no.model)) { # binary classification
      # tsp.list is a df here
      table <- tsp.list[1:K, ]
      kpairs[[as.character(K)]] <- rbind.data.frame(kpairs[[as.character(K)]], table)
    } else {
      # multi-class classification
      ifelse (no.model == 3, kpairs[[as.character(K)]] <- list(), kpairs[[as.character(K)]] <- data.frame())
      for (s in 1:length(tsp.list)) {
        name <- names(tsp.list)[s]
        scores <- tsp.list[[s]]
        if (no.model == 3) { # for model 3
          kpairs[[as.character(K)]][[name]] <- scores[1:K, ] %>% mutate(class = rep(name, K)) %>% select(-score)
        } else { # for model 1&2 merging
          table <- scores[1:K, ] %>% mutate(class = rep(name, K)) %>% select(-score)
          kpairs[[as.character(K)]] <- rbind.data.frame(kpairs[[as.character(K)]], table)
        }
      }
      # Delete repeated rows
      if (no.model != 3) {for (i in 1:length(kpairs)) {kpairs[[i]] <- delete.dup(kpairs[[i]])}}
    }
  }
  kpairs
}



