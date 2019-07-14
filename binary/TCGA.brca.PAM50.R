##############################################
## Call PAM50 Subtypes for TCGA 1095 patients
## Zhilin Pan
## July 2019
## Modified from the script of Zhuli 2016
##############################################

# Here I use the expression values which are median ratio normalized and variance stablizing transformed, 
# rather than TPM, since TPM is not cross-sampled normalized and thus not appropriate.


library(genefu)
library(genefilter)
data(pam50.robust)
data(pam50.scale)
data(pam50)

# Load data
tcga.counts <- get(load("~/ktsp/data/TCGA_vstFC_forPAM50.Rdata"))
clin <- get(load("~/ktsp/new.data/tcgaClinical_10_14_15_Processed.rda"))
rm(vst.counts, tcga_clinical_10_14_15) 

# Select samples
clin2 <- clin[match(colnames(tcga.counts), rownames(clin)), ]
ER <- clin2[ ,"er_status_by_ihc"]
ILC <- clin2[ ,"Histology_textCode"]
ILC[is.na(ILC)] <- "NA"
table(as.character(ER))
table(as.character(ER), as.character(ILC))
dim(clin2) #1095  130
dim(tcga.counts) #23368  1095


# Match gene names
gene50 <- rownames(pam50.robust$centroids.map)
length(intersect(rownames(tcga.counts), gene50)) #47
rownames(tcga.counts) <- sub("NUF2", "CDCA1", rownames(tcga.counts))
rownames(tcga.counts) <- sub("NDC80", "KNTC2", rownames(tcga.counts))
rownames(tcga.counts) <- sub("ORC6", "ORC6L", rownames(tcga.counts))
length(intersect(rownames(tcga.counts), gene50)) #50


# Subsample to make ER status balanced
cnt.erneg <- tcga.counts[match(gene50, rownames(tcga.counts)), which(ER=="Negative")]
dim(cnt.erneg)  #237
cnt.erpos <- tcga.counts[match(gene50, rownames(tcga.counts)), which(ER=="Positive")]
dim(cnt.erpos)  #808 (should sample 351 59.7%)


# Subsample to ensure same percentage of ERpos
set.seed(2013063)
subsample.erpos <- cnt.erpos[ ,sample(seq(1, ncol(cnt.erpos)), 237, replace=F)]
subsample <- cbind(cnt.erneg, subsample.erpos)
dim(subsample) #50 474


# Adjust to median50 --> purpose?
median50 <- apply(subsample, 1, median)
adj.counts <- tcga.counts[match(gene50, rownames(tcga.counts)), ]
adj.counts <- adj.counts - median50

# Predict (genefu) for one time
pam50.call <- intrinsic.cluster.predict(sbt.model = pam50, 
                                        data = t(adj.counts), 
                                        annot = annot.nkis, 
                                        do.mapping = FALSE, 
                                        do.prediction.strength = FALSE, 
                                        verbose = TRUE)
table(pam50.call$subtype)
table(pam50.call$subtype[which(ILC=="ILC")])
sum(ILC=="ILC",na.rm=T)  # LumA 80.7%


# Predict (compare correlation) --> ??
pam50cent <- pam50$centroids
assign <- rep(NA, ncol(adj.counts))
corre <- matrix(NA, ncol(adj.counts),5)
colnames(corre) <- c("Basal", "Her2", "LumA", "LumB", "Normal")
for(i in 1:ncol(adj.counts)){
  corre[i, ] <- cor(adj.counts[,i], pam50cent, method="spearman")
  assign[i] <- colnames(corre)[which.max(corre[i,])]
}
table(assign)
table(assign[which(ILC=="ILC")])

## compare
sum(assign != pam50.call$subtype)


# Use METABRIC Method (repeat 100 sampling) 
all.call <- matrix(NA, ncol(adj.counts), ncol = 100)
set.seed(2013063)
for(s in 1:100){
  # subsample to ensure same percentage of ERpos
  subsample.erpos <- cnt.erpos[ ,sample(seq(1, ncol(cnt.erpos)), 237, replace=F)]
  subsample <- cbind(cnt.erneg, subsample.erpos)
  median50.rep <- apply(subsample, 1, median)
  adj.counts.rep <- tcga.counts[match(gene50, rownames(tcga.counts)), ] - median50.rep
  pam50.call.sub <- intrinsic.cluster.predict(sbt.model = pam50, 
                                              data = t(adj.counts.rep), 
                                              annot = annot.nkis, 
                                              do.mapping = FALSE, 
                                              do.prediction.strength = FALSE, 
                                              verbose = TRUE)
  all.call[ ,s] <- pam50.call.sub$subtype 
  print(s)
}

final.assign.mul <- sapply(1:ncol(adj.counts), function(x) names(table(all.call[x, ])))
leng <- sapply(1:ncol(adj.counts), function(x) length(final.assign.mul[[x]]))
table(leng)

final.assign <- sapply(1:ncol(adj.counts), function(x) names(table(all.call[x, ]))[which.max(as.numeric(table(all.call[x, ])))])
names(final.assign) <- colnames(adj.counts)
table(final.assign)
table(final.assign[which(ILC=="ILC")])

# Compare results from 1 & 100 times
table(final.assign, pam50.call$subtype)
table(final.assign[which(ILC=="ILC")], pam50.call$subtype[which(ILC=="ILC")])


## Compare with the PAM50 assignments using TPM data
tpm.res <- get(load("~/ktsp/new.data/TCGA_PAM50_histology_TumorPurity.rda"))
table(final.assign, tpm.res$final_assign)

brca.clin <- cbind.data.frame(data.frame(final_assign=final.assign),
             clin2[ ,c('gender', 'race', 'ajcc_pathologic_tumor_stage', 'age_at_diagnosis',
                       'er_status_by_ihc', 'pr_status_by_ihc', 'her2_status_by_ihc', 'Histology_textCode')])
colnames(brca.clin) <- c('final_assign', 'gender', 'race', 'ajcc_tumor_stage', 'age_at_diagnosis',
                         'ER_status', 'PR_status', 'Her2_status', 'histology')

brca.clin$ajcc_tumor_stage <- gsub('^Stage ', '', brca.clin$ajcc_tumor_stage)
brca.clin$ajcc_tumor_stage <- gsub('[ABC]$', '', brca.clin$ajcc_tumor_stage)
brca.clin$ajcc_tumor_stage <- factor(brca.clin$ajcc_tumor_stage)

save(brca.clin, file='~/ktsp/data/BRCA.clin.Rdata')

