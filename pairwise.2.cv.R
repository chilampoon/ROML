suppressPackageStartupMessages({
  library(randomForest)
  library(dplyr)
  library(caret)
  library(pROC)
  library(ggpubr)
  library(reshape2)
  library(nloptr)
  library(stringr)
})

main.dir <- "~/GoogleDrive/pjs/kTSP"
data.dir <- file.path(main.dir, "data")
save.dir <- file.path(data.dir, "multi.results")
source(file.path(main.dir, "code/utils.R"))

# Load processed expression data
load(file.path(data.dir, "BRCA.expr.Rdata"))
load(file.path(data.dir, "BRCA.clinic.Rdata"))
load(file.path(data.dir, "MB.expr.Rdata"))
load(file.path(data.dir, "MB.clinic.Rdata"))


singleKpairs <- function(Ks, scores.dir) {
  kpairs <- list()
  all.rda <- list.files(path=scores.dir, pattern="*.scores.Rdata", full.names=T, recursive=F)
  
  # To remove sth. if exists (avoid memory boooomb)
  ifrm <- function(x, env = globalenv()) {
    if(exists(x, envir = env)) {
      rm(list = x, envir = env)
    }
  }
  
  
  for (rda in all.rda) {
    name <- gsub(".scores.Rdata$", "", rda)
    name <- strsplit(name, "/\\s*(?=[^/]+$)", perl=T)[[1]][2]
    kpairs[[name]] <- list()
    scores <- get(load(rda))
    scores$genes <- rownames(scores)
    try({
      for (K in Ks) {
        table <- scores[1:K, ]
        table <- table %>% mutate(geneX = gsub("^(.*?),.*", "\\1", table$genes), 
                                  geneY = gsub("^(.*?),", "", table$genes),
                                  class = rep(name, K)) %>% select(-score, -genes)
        
        kpairs[[name]][[as.character(K)]] <- table
      }
      
      ifrm("tmp") # pairwise objects names
      ifrm("scores") 
      ifrm(list=ls(pattern=".scores$")) # one vs other objects names
    })
  }
  kpairs
}

object.rf.cv <- function(nfold, Ks, labels, dataframe, kpairs) {
  # Create folds
  folds <- createFolds(labels, k = nfold)
  results.list <- list()
  
  # Outer loop for cross-validation
  # Intermediate loop for different Ks
  # Inner loop for different classes
  for (i in 1:5) {
    results.list[[i]] <- list()
    # 1/5 of the dataset as validation, remaining data for training
    cv.train <- dataframe[-folds[[i]], ]
    cv.train.labels <- labels[-folds[[i]]]
    cv.validate <- dataframe[folds[[i]], ]
    cv.validate.labels <- labels[folds[[i]]]
    
    for (K in Ks) {
      results.list[[i]][[as.character(K)]] <- list()
      grp.names <- names(kpairs)
      
      # Get RF models from each group pair
      prs.table <- data.frame()
      for (grp in grp.names) {
        group1 <- strsplit(grp, "_")[[1]][1]
        group2 <- strsplit(grp, "_")[[1]][2]
        ktsp <- kpairs[[grp]][[as.character(K)]]
        
        # Create binary matrices for gene ranks
        bi.cv.train <- cv.train[which(cv.train.labels == group1 | cv.train.labels == group2),]
        bi.cv.train.labels <- factor(cv.train.labels[cv.train.labels == group1 | cv.train.labels == group2])
        bi.cv.train <- bi.cv.train[, ktsp$geneX] - bi.cv.train[, ktsp$geneY]
        bi.cv.train <- as.data.frame(ifelse(bi.cv.train > 0, "1", "0"))
        colnames(bi.cv.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
        bi.cv.train <- cbind.data.frame(bi.cv.train, type=bi.cv.train.labels)
        
        bi.cv.validate <- cv.validate[, ktsp$geneX] - cv.validate[, ktsp$geneY]
        bi.cv.validate <- as.data.frame(ifelse(bi.cv.validate > 0, "1", "0"))
        colnames(bi.cv.validate) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
        bi.cv.validate <- cbind.data.frame(bi.cv.validate, type=cv.validate.labels)
        for (p in names(bi.cv.validate)) {
          if (p != "type") {
            levels(bi.cv.validate[[p]]) <- c("0", "1") -> levels(bi.cv.train[[p]])
          }
        } # Ensure the levels of predictors are same in train & validate
        
        # Build random forest for pair-wise groups
        rf.model <- randomForest(type~., data = bi.cv.train, ntree = 500, replace = F) # tuned later
        #predictionLabels <- predict(rf.model, newdata=bi.cv.validate, type='response')
        predictionPrs <- as.data.frame(predict(rf.model, newdata=bi.cv.validate, type='prob'))
        colnames(predictionPrs) <- paste0(colnames(predictionPrs), paste0(".", grp))
        
        # Combine & merge probability tables
        if (nrow(prs.table) == 0) {
          prs.table <- rbind.data.frame(prs.table, predictionPrs)
        } else {
          prs.table <- cbind.data.frame(prs.table, predictionPrs)
        }
      }
      
      # Optimize final probabilities using objective function
      opti.RFs <- function(prs.table, classes, pw.classes) {
        ################ Ensure the concordance of orders of classes and pw.classes ########################
        # Define an order for dataframe and x0
        pw.classes <- c(sprintf("%s.%s", unlist(lapply(strsplit(pw.classes, "_"), '[[', 1)), pw.classes),
                        sprintf("%s.%s", unlist(lapply(strsplit(pw.classes, "_"), '[[', 2)), pw.classes))
        prs.table <- prs.table[, pw.classes]
        
        lb <- rep(0, length(classes)) # lower bound
        ub <- rep(1, length(classes)) # upper bound
        x0 <- setNames(rep(1/length(classes), length(classes)), classes) # vector with starting values for the optimization
        
        # obj_f <- function(x, pw.p) {
        # def.var <- data.frame(class = c("Normal", "LumA", "LumB", "Her2", "Basal"),
        #                       numVar = seq(1, length(x)))
        # min.terms <- setNames(rep(0, length(colnames(pw.p))), colnames((pw.p)))
        # for (col in colnames(pw.p)) {
        #   class <- strsplit(col, "\\.")[[1]][1]
        #   pw.grp <- strsplit(col, "\\.")[[1]][2]
        #   pw.class <- strsplit(pw.grp, "_")[[1]][strsplit(pw.grp, "_")[[1]] != class]
        #   # Formula
        #   class <- def.var[def.var$class==class,]$numVar
        #   pw.class <- def.var[def.var$class==pw.class,]$numVar
        #   min.terms[col] <- (pw.p[, col] - (x0[class] / (x0[class] + x0[pw.class])))^2
        # } 
        # 
        # sum(min.terms)
        # }
        
        # Objective function
        obj_f <- function(x, pw.p) {
          # x0 Normal   LumA   LumB   Her2  Basal 
          #      0.2    0.2    0.2    0.2    0.2 
          
          (pw.p[,1] - x[4]/(x[4]+x[5]))^2 + (pw.p[,2] - x[2]/(x[2]+x[5]))^2 + (pw.p[,3] - x[2]/(x[2]+x[4]))^2 + 
            (pw.p[,4] - x[2]/(x[2]+x[3]))^2 + (pw.p[,5] - x[3]/(x[3]+x[5]))^2 + (pw.p[,6] - x[3]/(x[3]+x[4]))^2 + 
            (pw.p[,7] - x[1]/(x[1]+x[5]))^2 + (pw.p[,8] - x[1]/(x[1]+x[4]))^2 + (pw.p[,9] - x[1]/(x[1]+x[2]))^2 + 
            (pw.p[,10] - x[1]/(x[1]+x[3]))^2 + (pw.p[,11] - x[5]/(x[4]+x[5]))^2 + (pw.p[,12] - x[5]/(x[2]+x[5]))^2 + 
            (pw.p[,13] - x[4]/(x[2]+x[4]))^2 + (pw.p[,14] - x[3]/(x[2]+x[3]))^2 + (pw.p[,15] - x[5]/(x[5]+x[5]))^2 + 
            (pw.p[,16] - x[4]/(x[3]+x[4]))^2 + (pw.p[,17] - x[5]/(x[1]+x[5]))^2 + (pw.p[,18] - x[4]/(x[1]+x[4]))^2 + 
            (pw.p[,19] - x[2]/(x[1]+x[2]))^2 + (pw.p[,20] - x[3]/(x[1]+x[3]))^2
        }
        
        # Rewrite the equality as 1-(P1+P2+P3+P4+P5)=0
        eq_c <- function(x, pw.p) 1 - sum(x) 
      
        # Rewrite the inequality as P1*P2*P3*P4*P5 - 1 <= 0
        ineq_c <- function(x, pw.p) prod(x) - 1 
        
        # Summarize final probabilities
        final.prs <- setNames(data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors=F), classes)
        for (r in 1:nrow(prs.table)) {
          opt <- nloptr(x0=x0, eval_f=obj_f, lb=lb, ub=ub, eval_g_ineq=ineq_c, eval_g_eq=eq_c, 
                        opts=list(algorithm="NLOPT_GN_ISRES", "maxeval" = 100000), pw.p=prs.table[r, ])
          final.prs[nrow(final.prs)+1, ] <- opt$solution
        }
        rownames(final.prs) <- rownames(prs.table)
        final.prs
      }
      
      classes <- unique(as.character(labels))
      pw.classes <- names(sin.kpairs)
      final.prs <- opti.RFs(prs.table, classes, pw.classes)
      final.prs$predict.grp <- colnames(final.prs)[apply(final.prs, 1, which.max)] -> predictionLabels
      
      # Save performances
      cm <- confusionMatrix(data=as.factor(predictionLabels), reference=as.factor(cv.validate.labels))
      sensitivity <- as.numeric(cm$byClass[, 1])
      specificity <- as.numeric(cm$byClass[, 2])
      youden <- sensitivity + specificity - 1
      names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
      results.list[[i]][[as.character(K)]][['trueLabels']] <- cv.validate.labels
      results.list[[i]][[as.character(K)]][['predictionLabels']] <- predictionLabels
      results.list[[i]][[as.character(K)]][['predictionPrs']] <- final.prob
      results.list[[i]][[as.character(K)]][['ConfusionMatrix']] <- cm
      results.list[[i]][[as.character(K)]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
      results.list[[i]][[as.character(K)]][['Youden']] <- youden
    }
  
  }
  
}








########### MB as training; TCGA as testing
sub.genes <- lowMeans.filter(mb.expr, 25)
mb.train <- mb.expr[sub.genes, ]
brca.test <- brca.expr[sub.genes, ]


RF.lists <- list()
mb.scores.dir <- file.path(save.dir, "metabric/pairwise.scores")
K <- 20
labels <- factor(mb.clinic$NOT_IN_OSLOVAL_Pam50Subtype)
data.train <- as.data.frame(t(mb.train)) # Row - sample; Col - gene
sin.kpairs <- singleKpairs(K, mb.scores.dir)











Ks <- seq(10, 200, by = 10)
labels <- as.factor(brca.clinic$final_assign)
data.train <- t(brca.train) # X - sample; Y - gene
kpairs <- mergePairs(Ks, tcga.scores.dir)
cv.results <- tsp.rf.cv(5, Ks, labels, data.train, kpairs)


