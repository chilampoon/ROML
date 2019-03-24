#!/usr/bin/env Rscript

library(dplyr)
library(switchBox)

main.dir <- "~/GoogleDrive/pjs/kTSP"
data.dir <- file.path(main.dir, "data")
save.dir <- file.path(data.dir, "multi.results")
source(file.path(main.dir, "code/utils.R"))

# Load processed data
load(file.path(data.dir, "BRCA.expr.Rdata"))
load(file.path(data.dir, "BRCA.clinic.Rdata"))
load(file.path(data.dir, "MB.expr.Rdata"))
load(file.path(data.dir, "MB.clinic.Rdata"))

####### BRCA scores
sub.genes <- lowMeans.filter(brca.expr, 25)
brca.train <- brca.expr[sub.genes, ]
mb.test <- mb.expr[sub.genes, ]

# Subset BRCA into 5 expression sets according to subtypes
grp.names <- as.character(unique(brca.clinic$final_assign))
brca.train <- cbind.data.frame(t(brca.train), type=brca.clinic$final_assign)
sub.brca.expr <- list()
for (grp in grp.names) {
  sub.brca.expr[[grp]] <- brca.train[brca.train$type == grp, ]
}



# Get the tsp table for each paired group: Her2 vs LumA; Her2 vs LumB....
for (g in 1:length(grp.names)) {
  if (g != length(grp.names)) {
    try({
      for (s in (g+1):length(grp.names)) {
        grp.data <- rbind.data.frame(sub.brca.expr[[grp.names[g]]], sub.brca.expr[[grp.names[s]]])
        grp.vec <- factor(grp.data$type)
        grp.data <- t(as.matrix(grp.data[, which(colnames(grp.data) != "type")]))
        tmp <- SWAP.CalculateScores(grp.data, grp.vec, FilterFunc = NULL)
        tmp <- data.frame(score = sort(tmp$score, decreasing = T))
        #tmp$geneX <- gsub("^(.*?),.*", "\\1", rownames(tmp))
        #tmp$geneY <- gsub("^(.*?),", "", rownames(tmp))
        save(tmp, file=file.path(save.dir, paste0("tcga/",grp.names[g], "_", grp.names[s], ".scores.Rdata")))
        print(paste0("Finish ", paste0(grp.names[g], "_", grp.names[s]), "!!"))
      }
    })
  }
}

#save(tcga.pw.tables, file=file.path(save.dir, "tcga.50p.pw.tables.Rdata"))
## The list is too large to load :<


############ MetaBric scores

sub.genes2 <- lowMeans.filter(mb.expr, 25)
mb.train <- mb.expr[sub.genes2, ]

# Subset MB into 5 expression sets according to subtypes
grp.names <- as.character(unique(mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype))
mb.train <- cbind.data.frame(t(mb.train), type=mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype)
sub.mb.expr <- list()
for (grp in grp.names) {
  sub.mb.expr[[grp]] <- mb.train[mb.train$type == grp, ]
}


# Get the tsp table for each paired group: Her2 vs LumA; Her2 vs LumB....
for (g in 1:length(grp.names)) {
  if (g != length(grp.names)) {
    try({
      for (s in (g+1):length(grp.names)) {
        grp.data <- rbind.data.frame(sub.mb.expr[[grp.names[g]]], sub.mb.expr[[grp.names[s]]])
        grp.vec <- factor(grp.data$type)
        grp.data <- t(as.matrix(grp.data[, which(colnames(grp.data) != "type")]))
        tmp <- SWAP.CalculateScores(grp.data, grp.vec, FilterFunc = NULL)
        tmp <- data.frame(score = sort(tmp$score, decreasing = T))
        #tmp$geneX <- gsub("^(.*?),.*", "\\1", rownames(tmp))
        #tmp$geneY <- gsub("^(.*?),", "", rownames(tmp))
        save(tmp, file=file.path(save.dir, paste0("metabric/",grp.names[g], "_", grp.names[s], ".scores.Rdata")))
        print(paste0("Finish MB ", paste0(grp.names[g], "_", grp.names[s]), "!!"))
      }
    })
  }
}

