#!/usr/bin/env Rscript

library(switchBox)
library(dplyr)
source("~/ktsp/script/utils.R")

main.dir <- "~/ktsp"
data.dir <- file.path(main.dir, "data")
save.dir <- file.path(main.dir, "multi.results")

# Load processed data
load(file.path(data.dir, "BRCA.expr.Rdata"))
load(file.path(data.dir, "BRCA.clinic.Rdata"))
load(file.path(data.dir, "MB.expr.Rdata"))
load(file.path(data.dir, "MB.clinic.Rdata"))

#### Calculate TSP scores (modified from Kelly's script)
# calTSPscores <- function(data, grp.vec) {
#   t.data <- t(data)
#   n <- dim(t.data)[1] # number of samples
#   m <- dim(t.data)[2] # number of genes
#   xx <- as.matrix(cbind(rep(1, n), grp.vec))
#   hatxx <- solve(t(xx) %*% xx) %*% t(xx) # '%*%': matrix multiplication; solve(): return the inverse of the matrix
#   
#   # Calculate paired-gene distance
#   print("  Calculate paired-gene distance..")
#   registerDoParallel()
#   tmp <- foreach(j = c(1:(m-1))) %dopar% {
#     yy <- t.data[ ,-c(1:j)] < t.data[ ,j] # n x m matrix #compares every gene to each other non-repeatly (pairs)
#     t(t((hatxx %*% yy)[2, ]))
#   }
#   out <- do.call(rbind, tmp)
#   score <- rowSums(out)
#   names(score) <- 1:length(score)
#   order.score <- score[order(-abs(score))] # note that using absolute values here!!
#   tmp.thread <- as.numeric(names(order.score))
#   
#   # Assign names of gene pairs
#   ind <- foreach(k = 1:(m-1)) %dopar% {
#     t.data  <- t(data)
#     #cbind(rep(k, length(k:m)-1), (k+1):m, rep(colnames(t.data)[k], length(k:m)-1), colnames(t.data)[(k+1):m])
#     cbind(rep(colnames(t.data)[k], length(k:m)-1), colnames(t.data)[(k+1):m])
#   }
#   df <- do.call(rbind, ind)
#   df <- df[tmp.thread, ]
#   df <- cbind(df, order.score)
#   colnames(df) <- c("geneX", "geneY", "score")
#   df <- as.data.frame(df)
#   df
# }


## TCGA-BRCA scores
sub.brca <- lowMeans.filter(brca.expr, 25)
brca.train <- brca.expr[sub.brca, ]


# Get the tsp table for each outcome
all(colnames(brca.train) == rownames(brca.clinic)) # T
grp <- brca.clinic$final_assign
grp.names <- levels(grp)


print("Start calculating BRCA tsp scores...")
for (group in grp.names) {
  try({
    grp.vec <- rep(0, length(grp))
    grp.vec <- factor(replace(grp.vec, which(grp == group), 1)) # one v.s. all other
    tmp <- SWAP.CalculateScores(as.matrix(brca.train), grp.vec, FilterFunc = NULL)
    tmp <- data.frame(score = sort(tmp$score, decreasing = T))
    tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
    tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
    rownames(tmp) <- NULL
    save(tmp, file = file.path(save.dir, paste0("tcga/one.other.scores/", group, ".score.Rdata")))
    print(paste0("Finish BRCA.", group, "!!"))
  })
}



## MetaBric scores
sub.mb <- lowMeans.filter(mb.expr, 25)
mb.train <- mb.expr[sub.mb, ]


# Get the tsp table for each outcome
mb.grp <- mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype
mb.names <- grp.names # concordance as before

for (group in mb.names) {
  try({
    grp.vec <- rep(0, length(mb.grp))
    grp.vec <- factor(replace(grp.vec, which(mb.grp == group), 1)) # one v.s. all other
    tmp <- SWAP.CalculateScores(as.matrix(mb.train), grp.vec, FilterFunc = NULL)
    tmp <- data.frame(score = sort(tmp$score, decreasing = T))
    tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
    tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
    rownames(tmp) <- NULL
    save(tmp, file = file.path(save.dir, paste0("metabric/one.other.scores/", group, ".score.Rdata")))
    print(paste0("Finisph MB.", group, "!!"))
  })
}

