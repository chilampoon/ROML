##############################################
## Call PAM50 Subtypes for TCGA 1095 patients
## Li Zhu
## 01/23/2016
##############################################
rm(list = ls())
#source("http://bioconductor.org/biocLite.R")
#biocLite("genefu")
library("genefu")
library("genefilter")
data(pam50.robust)
data(pam50.scale)
data(pam50)

setwd("/net/wong05/home/liz86/Steffi/Kevin_IDC_ILC_DE/Data")
load("brcaTCGAtpm_log2_uniquePrimaries.Rda")
data<-brcaTCGAtpm_log2_primaries_unique
load("tcgaClinical_10_14_15_Processed.Rda")
clin<-tcga_clinical_10_14_15

## select samples
clin_match<-clin[match(colnames(data),rownames(clin)),]
ER<-clin_match[,"er_status_by_ihc"]
ILC<-clin_match[,"Histology_textCode"]
ILC[is.na(ILC)]<-"NA"
table(as.character(ER))
clin2<-clin_match
data2<-data
table(as.character(ER),as.character(ILC))
dim(clin2) #1095  130
dim(data2) #23368  1095

## match gene names
gene50<-rownames(pam50.robust$centroids.map)
length(intersect(rownames(data2),gene50)) #47
rownames(data2) <- sub("NUF2", "CDCA1", rownames(data2))
rownames(data2) <- sub("NDC80", "KNTC2", rownames(data2))
rownames(data2) <- sub("ORC6", "ORC6L", rownames(data2))
length(intersect(rownames(data2),gene50)) #50

## subsample to make ER status balanced
data_erneg<-data2[match(gene50,rownames(data2)),which(ER=="Negative")]
dim(data_erneg)  #237
data_erpos<-data2[match(gene50,rownames(data2)),which(ER=="Positive")]
dim(data_erpos)  #808 (should sample 351 59.7%)

# subsample to ensure same percentage of ERpos
set.seed(2013063)
subsample_erpos<-data_erpos[,sample(seq(1,ncol(data_erpos)),237,replace=F)]
subsample<-cbind(data_erneg,subsample_erpos)
dim(subsample) #50 474

# subsample to ensure same number of ERpos and ERneg (turns out to be worse)
if(F){
  set.seed(2013063)
  subsample_erneg<-data_erpos[,sample(seq(1,ncol(data_erpos)),77,replace=F)]
  subsample_erpos<-data_erpos[,sample(seq(1,ncol(data_erpos)),114,replace=F)]
  subsample<-cbind(subsample_erneg,subsample_erpos)
  dim(subsample) #50 191
}

median50<-apply(subsample,1,median)

## adjust to median50
data3<-data2[match(gene50,rownames(data2)),]
data3_c<-data3-median50

## predict (genefu)
pam50_call<-intrinsic.cluster.predict(sbt.model=pam50, 
                                      data=t(data3_c), 
                                      annot=annot.nkis, 
                                      do.mapping=FALSE, 
                                      do.prediction.strength=FALSE, 
                                      verbose=TRUE)
table(pam50_call$subtype)
table(pam50_call$subtype[which(ILC=="ILC")])
length(pam50_call$subtype)
sum(ILC=="ILC",na.rm=T)  # LumA 80.7%
save(pam50_call,file="pam50_call_TCGA_subsample_1095.RData")


## predict (compare correlation)
pam50cent<-pam50$centroids
assign<-rep(NA,ncol(data3_c))
corre<-matrix(NA,ncol(data3_c),5)
colnames(corre)<-c("Basal","Her2","LumA","LumB","Normal")
for(i in 1:ncol(data3_c)){
  corre[i,]<-cor(data3_c[,i],pam50cent,method="spearman")
  assign[i]<-colnames(corre)[which.max(corre[i,])]
}
table(assign)
table(assign[which(ILC=="ILC")])

## compare
sum(assign!=pam50_call$subtype)

### Use METABRIC Method (repeat 100 sampling) ###
all_call<-matrix(NA,ncol(data3),100)
set.seed(2013063)
for(s in 1:100){
  # subsample to ensure same percentage of ERpos
  subsample_erpos<-data_erpos[,sample(seq(1,ncol(data_erpos)),237,replace=F)]
  subsample<-cbind(data_erneg,subsample_erpos)
  median50_rep<-apply(subsample,1,median)
  data3_c_rep<-data2[match(gene50,rownames(data2)),]-median50_rep
  pam50_call_sub<-intrinsic.cluster.predict(sbt.model=pam50, 
                                            data=t(data3_c_rep), 
                                            annot=annot.nkis, 
                                            do.mapping=FALSE, 
                                            do.prediction.strength=FALSE, 
                                            verbose=TRUE)
  all_call[,s]<-pam50_call_sub$subtype 
  print(s)
}
final_assign_mul<-sapply(1:ncol(data3),function(x)names(table(all_call[x,])))
leng<-sapply(1:ncol(data3),function(x)length(final_assign_mul[[x]]))
table(leng)

final_assign<-sapply(1:ncol(data3),function(x)names(table(all_call[x,]))[which.max(as.numeric(table(all_call[x,])))])
names(final_assign)<-colnames(data3)
table(final_assign)
table(final_assign[which(ILC=="ILC")])

table(final_assign,pam50_call$subtype)
table(final_assign[which(ILC=="ILC")],pam50_call$subtype[which(ILC=="ILC")])

save(final_assign,file="PAM50_TCGA_useMetabricMethod_1095.RData")

table(final_assign[which(ILC2 %in% c("ILC","IDC"))],ILC2[which(ILC2 %in% c("ILC","IDC"))])

### compare with TCGA ILC paper
paper<-read.csv("Comprehensive_ILC_fig1.csv",skip=2)

## compare subjects
length(intersect(paper[,1],rownames(clin)))  #817
clin_match_paper<-clin[match(paper[,1],rownames(clin)),]

## compare ILC labels
ILC_label_paper<-paper[,"Final.Pathology"]
ILC_label_our<-clin_match_paper[,"Histology_textCode"]
ILC_label_our[is.na(ILC_label_our)]<-"NA"
table(as.character(ILC_label_our),as.character(ILC_label_paper))

## compare ER status
ER_label_paper<-paper[,"ER.IHC"]
ER_label_our<-clin_match_paper[,"er_status_by_ihc"]
table(as.character(ER_label_our),as.character(ER_label_paper))

## compare PAM50
ourpam50<-final_assign[match(paper[,1],colnames(data3))]
table(ourpam50,paper[,3])

paper_ILC<-paper[which(paper[,2]=="ILC"),]
ourILC<-final_assign[match(paper_ILC[,1],colnames(data3))]

length(intersect(paper_ILC[,1],colnames(data3)))  #127
ourpam50_ILC<-final_assign[match(paper_ILC[,1],colnames(data3))]
table(ourpam50_ILC,paper_ILC[,3])

ourpam50_ILC_nc<-pam50_call_not_centered$subtype[match(paper_ILC[,1],colnames(data3))]
table(ourpam50_ILC_nc,paper_ILC[,3])

# metabric method
table(final_assign,pam50_call$subtype)
ourpam50<-final_assign[match(paper[,1],colnames(data3))]
table(ourpam50,paper[,3])
ourpam50_ILC<-final_assign[match(paper_ILC[,1],colnames(data3))]
table(ourpam50_ILC,paper_ILC[,3])


