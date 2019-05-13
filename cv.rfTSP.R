# Functions to do the cross-validation
suppressPackageStartupMessages({
  library(randomForest)
  library(nloptr)
  library(caret)
})



tsp.rf.cv <- function(nfold, Ks, labels, dataframe, grp.names, kpairs, no.model, baseline=FALSE) {
  # Create folds
  folds <- createFolds(labels, k=nfold)
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
    
    if (baseline == TRUE) { 
      # Baseline model
      cv.train <- cbind.data.frame(cv.train, type=cv.train.labels)
      cv.validate <- cbind.data.frame(cv.validate, type=cv.validate.labels)
      
      # Build the model
      rf.model <- randomForest(type~., data=cv.train, ntree=500, replace=F)
      predictionLabels <- predict(rf.model, newdata=cv.validate, type='response')
      predictionScores <- predict(rf.model, newdata=cv.validate, type='prob')
      
      # Save performances
      cm <- confusionMatrix(data=as.factor(predictionLabels), reference=as.factor(cv.validate.labels))
      sensitivity <- as.numeric(cm$byClass[, 1])
      specificity <- as.numeric(cm$byClass[, 2])
      youden <- sensitivity + specificity - 1
      names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
      results.list[[i]][['trueLabels']] <- cv.validate.labels
      results.list[[i]][['predictionLabels']] <- predictionLabels
      results.list[[i]][['predictionScores']] <- predictionScores
      results.list[[i]][['ConfusionMatrix']] <- cm
      results.list[[i]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
      results.list[[i]][['Youden']] <- youden
      results.list[[i]][['RFmodels']] <- rf.model
      #print(paste0("Fold ", i, " ACC:", round(results.list[[i]][['ACC']], 5)))
    } else {
      for (K in Ks) {
        results.list[[i]][[as.character(K)]] <- list()
        #results.list[[i]][[as.character(K)]][['kpairs']] <- all.tsp
        
        # For model 1 & 2
        if (no.model != 3) { 
          all.tsp <- kpairs[[as.character(no.model)]][[as.character(K)]]
          # Create binary matrices for gene ranks
          bi.cv.train <- cv.train[, all.tsp$geneX] - cv.train[, all.tsp$geneY]
          bi.cv.train <- as.data.frame(ifelse(bi.cv.train > 0, "1", "0"))
          colnames(bi.cv.train) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
          bi.cv.train <- cbind.data.frame(bi.cv.train, type=cv.train.labels)
          
          bi.cv.validate <- cv.validate[, all.tsp$geneX] - cv.validate[, all.tsp$geneY]
          bi.cv.validate <- as.data.frame(ifelse(bi.cv.validate > 0, "1", "0"))
          colnames(bi.cv.validate) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
          bi.cv.validate <- cbind.data.frame(bi.cv.validate, type=cv.validate.labels)
          for (p in names(bi.cv.validate)) { # Ensure the levels of predictors are same in train & validate
            if (p != "type") levels(bi.cv.validate[[p]]) <- c("0", "1") -> levels(bi.cv.train[[p]])
          }
          
          # Build random forest
          rf.model <- randomForest(type~., data = bi.cv.train, ntree = 500, replace = F) 
          results.list[[i]][[as.character(K)]][['RFmodels']] <- rf.model
          predictionLabels <- predict(rf.model, newdata=bi.cv.validate, type='response')
          predictionPrs <- predict(rf.model, newdata=bi.cv.validate, type='prob')
        } else { # For model 3
          # Get RF models from each group pair
          #print('  Building RF models from each group pair')
          results.list[[i]][[as.character(K)]][['RFmodels']] <- list()
          prs.table <- data.frame()
          for (g in 1:(length(grp.names)-1)) {
            for (r in (g+1):length(grp.names)) {
              group1 <- grp.names[g]
              group2 <- grp.names[r]
              ktsp <- kpairs[[as.character(no.model)]][[as.character(K)]][[sprintf('%s_%s', group1, group2)]]
              
              # Create binary matrices for gene ranks
              bi.cv.train <- cv.train[which(cv.train.labels == group1 | cv.train.labels == group2),]
              bi.cv.train.labels <- factor(cv.train.labels[cv.train.labels == group1 | cv.train.labels == group2])
              bi.cv.train <- bi.cv.train[, ktsp$geneX] - bi.cv.train[, ktsp$geneY]
              bi.cv.train <- as.data.frame(ifelse(bi.cv.train > 0, "1", "0"))
              colnames(bi.cv.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
              bi.cv.train <- cbind.data.frame(bi.cv.train, type=bi.cv.train.labels)
              
              bi.cv.validate <- cv.validate[ ,ktsp$geneX] - cv.validate[ ,ktsp$geneY]
              bi.cv.validate <- as.data.frame(ifelse(bi.cv.validate > 0, "1", "0"))
              colnames(bi.cv.validate) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
              bi.cv.validate <- cbind.data.frame(bi.cv.validate, type=cv.validate.labels)
              for (p in names(bi.cv.validate)) {
                if (p != "type") {
                  levels(bi.cv.validate[[p]]) <- c("0", "1") -> levels(bi.cv.train[[p]])
                }
              } 
              
              # Build random forest for pairwise groups
              rf.model <- randomForest(type~., data = bi.cv.train, ntree = 500, replace = F) # tuned later
              results.list[[i]][[as.character(K)]][['RFmodels']][[sprintf('%s_%s', group1, group2)]] <- rf.model
              predictionPrs <- as.data.frame(predict(rf.model, newdata=bi.cv.validate, type='prob'))
              predictionPrs <- predictionPrs[ ,c(group1, group2)]
              colnames(predictionPrs) <- paste0(colnames(predictionPrs), paste0(".", sprintf('%s_%s', group1, group2)))
              
              # Combine & merge probability tables
              if (nrow(prs.table) == 0) {
                prs.table <- rbind.data.frame(prs.table, predictionPrs)
              } else {
                prs.table <- cbind.data.frame(prs.table, predictionPrs)
              }
            }
          }
          #print('  Calculating final probs using our obj func.')
          predictionPrs <- opti.RFs(prs.table, grp.names)
          predictionPrs$predict.grp <- colnames(predictionPrs)[apply(predictionPrs, 1, which.max)] -> predictionLabels
        }
        
        # Store performances
        predictionLabels <- factor(predictionLabels, levels = levels(labels))
        cm <- confusionMatrix(data=as.factor(predictionLabels), reference=as.factor(cv.validate.labels))
        sensitivity <- as.numeric(cm$byClass[, 1])
        specificity <- as.numeric(cm$byClass[, 2])
        youden <- sensitivity + specificity - 1
        names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
        results.list[[i]][[as.character(K)]][['trueLabels']] <- cv.validate.labels
        results.list[[i]][[as.character(K)]][['predictionLabels']] <- predictionLabels
        results.list[[i]][[as.character(K)]][['predictionScores']] <- predictionPrs
        results.list[[i]][[as.character(K)]][['ConfusionMatrix']] <- cm
        results.list[[i]][[as.character(K)]][['Sens']] <- sensitivity
        results.list[[i]][[as.character(K)]][['Specs']] <- specificity
        results.list[[i]][[as.character(K)]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
        results.list[[i]][[as.character(K)]][['Youden']] <- youden
        if (no.model != 3) results.list[[i]][[as.character(K)]][['RFmodels']] <- rf.model
        #print(paste0("Fold ", i, " & K=", K,  " ACC:", round(results.list[[i]][[as.character(K)]][['ACC']], 5)))
      }
    }
  }
  results.list
}



# Optimize final probabilities using objective function for model 3
opti.RFs <- function(prs.table, grp.names) {
  ngrp <- length(grp.names)
  lb <- rep(0, ngrp) # lower bound
  ub <- rep(1, ngrp) # upper bound
  x0 <- setNames(rep(1/ngrp, ngrp), grp.names) # vector with starting values for the optimization
  
  # Objective function
  obj_f <- function(x, pw.p.mat) {
    # x0 Normal   LumA   LumB   Her2  Basal 
    #      0.2    0.2    0.2    0.2    0.2 
    fml <- 0
    for (g in 1:(length(x)-1)) {
      for (r in (g+1):length(x)) {
        fml <- fml + (pw.p.mat[g, r] - x[g]/(x[g] + x[r]))^2 + (pw.p.mat[r, g] - x[r]/(x[g] + x[r]))^2
      }
    }
    sum(fml)
  }
  
  # Rewrite the equality as 1 - (P1 + P2 + P3 + P4 + P5) = 0
  eq_c <- function(x, pw.p.mat) 1 - sum(x) 
  
  # Rewrite the inequality as P1*P2*P3*P4*P5 - 1 <= 0
  ineq_c <- function(x, pw.p.mat) prod(x) - 1
  
  # Function to transform a row in prs.table to a matrix
  rowToMat <- function(row, grp.names) {
    row.mat <- setNames(data.frame(matrix(ncol=ngrp, nrow=ngrp), stringsAsFactors=F), grp.names)
    rownames(row.mat) <- colnames(row.mat)
    for (g in 1:(ngrp-1)) {
      grp1 <- grp.names[g]
      for (r in (g+1):ngrp) {
        grp2 <- grp.names[r]
        grp <- sprintf('%s_%s', grp1, grp2)
        row.mat[g, r] <- row[, sprintf('%s.%s', grp1, grp)]
        row.mat[r, g] <- row[, sprintf('%s.%s', grp2, grp)]
      }
    }
    row.mat
  }
  
  
  # Summarize final probabilities
  final.prs <- setNames(data.frame(matrix(ncol = ngrp, nrow = 0), stringsAsFactors=F), grp.names)
  for (n in 1:nrow(prs.table)) {
    pw.p.mat <- rowToMat(prs.table[n, ], grp.names)
    opt <- nloptr(x0=x0, eval_f=obj_f, lb=lb, ub=ub, eval_g_ineq=ineq_c, eval_g_eq=eq_c, 
                  opts=list(algorithm="NLOPT_GN_ISRES", "maxeval" = 150000), pw.p.mat=pw.p.mat)
    final.prs[nrow(final.prs)+1, ] <- opt$solution
  }
  rownames(final.prs) <- rownames(prs.table)
  final.prs
}



# Get overall values of performances from CV results
getOverall <- function(results.list, Ks, baseline=FALSE) {
  overall <- list()
  
  if (baseline == TRUE) {
    all.true.labels <- c()
    all.predict.labels <- c()
    all.predict.scores <- c()
    for (f in 1:length(results.list)) {
      all.true.labels <- c(all.true.labels, as.character(results.list[[f]][['trueLabels']]))
      all.predict.labels <- c(all.predict.labels, as.character(results.list[[f]][['predictionLabels']]))
      all.predict.scores <- rbind(all.predict.scores, results.list[[f]][['predictionScores']])
    }
    
    all.cm <- confusionMatrix(data=as.factor(all.predict.labels), reference=as.factor(all.true.labels))
    overall[['ACC']] <- as.numeric(all.cm$overall['Accuracy'])
    overall[['Sensitivity']] <- as.numeric(all.cm$byClass[, 1]) -> all.sens
    overall[['Specificity']] <- as.numeric(all.cm$byClass[, 2]) -> all.spec
    all.youden <- all.sens + all.spec - 1
    names(all.youden) <- gsub('Class: ', '', rownames(all.cm$byClass))
    overall[['Youden']] <- all.youden
    overall[['predictionScores']] <- all.predict.scores
    overall[['cm']] <- all.cm
    overall[['trueLabels']] <- all.true.labels
    overall[['avg.sens']] <- sum(all.sens)/length(all.sens)
    overall[['avg.spec']] <- sum(all.spec)/length(all.spec)
    overall[['avg.yd']] <- sum(all.youden)/length(all.youden)
  } else {
    for (k in 1:length(Ks)) {
      K <- Ks[k]
      overall[[as.character(K)]] <- list()
      all.true.labels <- c()
      all.predict.labels <- c()
      all.predict.scores <- c()
      for (f in 1:length(results.list)) {
        all.true.labels <- c(all.true.labels, as.character(results.list[[f]][[as.character(K)]][['trueLabels']]))
        all.predict.labels <- c(all.predict.labels, as.character(results.list[[f]][[as.character(K)]][['predictionLabels']]))
        all.predict.scores <- rbind(all.predict.scores, results.list[[f]][[as.character(K)]][['predictionScores']])
      }
      
      all.cm <- confusionMatrix(data=as.factor(all.predict.labels), reference=as.factor(all.true.labels))
      overall[[as.character(K)]][['ACC']] <- as.numeric(all.cm$overall['Accuracy'])
      overall[[as.character(K)]][['Sensitivity']] <- as.numeric(all.cm$byClass[, 1]) -> all.sens
      overall[[as.character(K)]][['Specificity']] <- as.numeric(all.cm$byClass[, 2]) -> all.spec
      all.youden <- all.sens + all.spec - 1
      names(all.youden) <- gsub('Class: ', '', rownames(all.cm$byClass))
      overall[[as.character(K)]][['Youden']] <- all.youden
      overall[[as.character(K)]][['predictionScores']] <- all.predict.scores
      overall[[as.character(K)]][['cm']] <- all.cm
      overall[[as.character(K)]][['avg.sens']] <- sum(all.sens)/length(all.sens)
      overall[[as.character(K)]][['avg.spec']] <- sum(all.spec)/length(all.spec)
      overall[[as.character(K)]][['avg.yd']] <- sum(all.youden)/length(all.youden)
    }
  }
  overall
}
