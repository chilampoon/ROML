# Preprocessing TCGA & MetaBric data
# Data downloaded from ktsp box

main.dir <- "~/ktsp"
data.dir <- file.path(main.dir, "data")
save.dir <- file.path(main.dir, "data")

# TCGA 
brca.expr <- get(load(file.path(data.dir, "brcaTCGAtpm_log2_uniquePrimaries.Rda")))
brca.clinic <- get(load(file.path(data.dir, "TCGA_PAM50_histology_TumorPurity.rda")))


# MetaBric
mb.expr <- get(load(file.path(data.dir, "MetaBric_processed_Li.RData")))
mb.clinic <- get(load(file.path(data.dir, "Complete_METABRIC_Clinical_Features_Data (1).rbin")))
rm(brcaTCGAtpm_log2_primaries_unique, TCGA_PAM50, MetaBric, Complete_METABRIC_Clinical_Features_Data)


# Intersect with two datasets
dim(brca.expr); dim(mb.expr)
common.genes <- intersect(rownames(brca.expr), rownames(mb.expr))
brca.expr <- brca.expr[common.genes, ]
mb.expr <- mb.expr[common.genes, ]


# Filter out MB NC samples
table(brca.clinic$final_assign); table(mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype)
ex_sample <- rownames(mb.clinic[mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype == "NC", ])
mb.expr <- mb.expr[, which(!colnames(mb.expr) %in% ex_sample)]
mb.clinic <- mb.clinic[which(!rownames(mb.clinic) %in% ex_sample), ]
mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype <- factor(mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype)
dim(brca.expr); dim(mb.expr)


# Correct the gene names
rownames(brca.expr) <- gsub("-", "_", rownames(brca.expr))
rownames(mb.expr) <- gsub("-", "_", rownames(mb.expr))


# Delete normal-like tumoral samples (4/23 update)
brca.clinic <- brca.clinic[brca.clinic$final_assign != 'Normal', ]
brca.clinic$final_assign <- factor(brca.clinic$final_assign)
dim(brca.clinic) # 1060 5
brca.expr <- brca.expr[, match(rownames(brca.clinic), colnames(brca.expr))]
dim(brca.expr) # 18202 1060
all(rownames(brca.clinic) == colnames(brca.expr)) # TRUE

mb.clinic <- mb.clinic[mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype != 'Normal', ]
mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype <- factor(mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype)
dim(mb.clinic) # 1775 25
mb.expr <- mb.expr[, match(rownames(mb.clinic), colnames(mb.expr))]
dim(mb.expr) # 18202  1775
all(rownames(mb.clinic) == colnames(mb.expr)) # TRUE


# Save & delete
save(brca.expr, file=file.path(save.dir, "BRCA.expr.Rdata"))
save(mb.expr, file=file.path(save.dir, "MB.expr.Rdata"))
save(brca.clinic, file=file.path(save.dir, "BRCA.clinic.Rdata"))
save(mb.clinic, file=file.path(save.dir, "MB.clinic.Rdata"))

file.remove(file.path(data.dir, "brcaTCGAtpm_log2_uniquePrimaries.Rda"))
file.remove(file.path(data.dir, "TCGA_PAM50_histology_TumorPurity.rda"))
file.remove(file.path(data.dir, "MetaBric_processed_Li.RData"))
file.remove(file.path(data.dir, "Complete_METABRIC_Clinical_Features_Data (1).rbin"))
