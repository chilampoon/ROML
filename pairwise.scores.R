#!/usr/bin/env Rscript

library(dplyr)
library(switchBox)

main.dir <- "~/ktsp"
data.dir <- file.path(main.dir, "data")
save.dir <- file.path(main.dir, "multi.results")
source("~/ktsp/script/utils.R")


# Load processed data
load(file.path(data.dir, "BRCA.expr.Rdata"))
load(file.path(data.dir, "BRCA.clinic.Rdata"))
load(file.path(data.dir, "MB.expr.Rdata"))
load(file.path(data.dir, "MB.clinic.Rdata"))


## BRCA scores
# Subset data
sub.brca <- lowMeans.filter(brca.expr, 25)
brca.train <- brca.expr[sub.brca, ]

# Get the tsp table for each outcome
all(colnames(brca.train) == rownames(brca.clinic)) # T
grp <- brca.clinic$final_assign
grp.names <- levels(grp)

# Subset BRCA into 5 expression sets according to subtypes
brca.train <- cbind.data.frame(t(brca.train), type=grp)
sub.brca.expr <- list()
for (grp in grp.names) {
  sub.brca.expr[[grp]] <- brca.train[brca.train$type == grp, ]
}

# Get the tsp table for each paired group: Her2 vs LumA; Her2 vs LumB....
for (i in 1:length(grp.names)) {
  if (i != length(grp.names)) {
    for (j in (i+1):length(grp.names)) {
      grp.data <- rbind.data.frame(sub.brca.expr[[grp.names[i]]], sub.brca.expr[[grp.names[j]]])
      grp.vec <- factor(grp.data$type)
      grp.data <- t(as.matrix(grp.data[, which(colnames(grp.data) != "type")]))
      tmp <- SWAP.CalculateScores(grp.data, grp.vec, FilterFunc = NULL)
      tmp <- data.frame(score = sort(tmp$score, decreasing = T))
      tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
      tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
      save(tmp, file = file.path(save.dir, paste0("tcga/pairwise.scores/", grp.names[i], "_", grp.names[j], ".score.Rdata")))
      print(paste0("Finish TCGA ", paste0(grp.names[i], "_", grp.names[j]), "!!"))
    }
  }
}


## MetaBric scores
sub.mb <- lowMeans.filter(mb.expr, 25)
mb.train <- mb.expr[sub.mb, ]

# Get the tsp table for each outcome
mb.grp <- mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype
mb.names <- grp.names # concordance as before

for (i in 1:length(mb.names)) {
  if (i != length(mb.names)) {
    for (j in (i+1):length(mb.names)) {
      grp.data <- rbind.data.frame(sub.brca.expr[[mb.names[i]]], sub.brca.expr[[mb.names[j]]])
      grp.vec <- factor(grp.data$type)
      grp.data <- t(as.matrix(grp.data[, which(colnames(grp.data) != "type")]))
      tmp <- SWAP.CalculateScores(grp.data, grp.vec, FilterFunc = NULL)
      tmp <- data.frame(score = sort(tmp$score, decreasing = T))
      tmp$geneX <- gsub('^(.*?),.*', '\\1', rownames(tmp))
      tmp$geneY <- gsub('^(.*?),', '', rownames(tmp))
      save(tmp, file = file.path(save.dir, paste0("metabric/pairwise.scores/", grp.names[i], "_", grp.names[j], ".score.Rdata")))
      print(paste0("Finish MB ", paste0(grp.names[i], "_", grp.names[j]), "!!"))
    }
  }
}

