# Functions to do the cross-validation
suppressPackageStartupMessages({
  library(randomForest)
  library(caret)
  library(caTools)
  library(pROC)
})
source('~/ktsp/script/function/objective.fun.R')

multi.rf.cv <- function(dataframe, label, K, kpairs, grp.names, model.no, time=10) {
  # Create folds
  if (nrow(dataframe) != length(label)) stop('Matrix input should be a dataframe which is sample by gene.')
  folds <- createFolds(label, k=5)
  
  # to final benchmarking
  all.acc <- c()
  all.spec <- data.frame() -> all.sens -> all.youden
  all.cv <- list()
  
  # repeat 10 times
  for (t in 1:time) {
    #print(paste0('Time=', t))
    all.cv[[as.character(t)]] <- list() -> cv.result
  
    # Outer loop for cross-validation
    # Inner loop for different classes
    for (i in 1:5) {
      cv.result[[i]] <- list()
      # 1/5 of the dataset as validation, remaining data for training
      train <- dataframe[-folds[[i]], ]
      train.label <- label[-folds[[i]]]
      valid <- dataframe[folds[[i]], ]
      valid.label <- label[folds[[i]]]
      
      if (missing(K)) {
        # For baseline model
        cv.train <- cbind.data.frame(train, type=train.label)
        cv.valid <- cbind.data.frame(valid, type=valid.label)
        
        # Build the model
        rf.model <- randomForest(type~., data=cv.train, ntree=501, replace=F)
        predictionLabel <- predict(rf.model, newdata=cv.valid, type='response')
        predictionScore <- predict(rf.model, newdata=cv.valid, type='prob')
      } else {
        # For model 1 & 2
        if (model.no != 3) {
          all.tsp <- kpairs[[as.character(K)]]
          # Create binary matrices for gene ranks
          cv.train <- train[ ,all.tsp$geneX] - train[ ,all.tsp$geneY]
          cv.train <- as.data.frame(ifelse(cv.train > 0, "1", "0"))
          colnames(cv.train) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
          cv.train <- cbind.data.frame(cv.train, type=train.label)
          
          cv.valid <- valid[ ,all.tsp$geneX] - valid[ ,all.tsp$geneY]
          cv.valid <- as.data.frame(ifelse(cv.valid > 0, "1", "0"))
          colnames(cv.valid) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
          cv.valid <- cbind.data.frame(cv.valid, type=valid.label)
          for (p in names(cv.valid)) { # Ensure the levels of predictors are same in train & validate
            if (p != "type") levels(cv.valid[[p]]) <- c("0", "1") -> levels(cv.valid[[p]])
          }
          # Build the model
          rf.model <- randomForest(type~., data=cv.train, ntree=501, replace=F)
          predictionLabel <- predict(rf.model, newdata=cv.valid, type='response')
          predictionScore <- predict(rf.model, newdata=cv.valid, type='prob')
        } else { 
          # For model 3
          # Get RF models from each group pair
          prs.table <- data.frame()
          for (g in 1:(length(grp.names)-1)) {
            for (r in (g+1):length(grp.names)) {
              group1 <- grp.names[g]
              group2 <- grp.names[r]
              ktsp <- kpairs[[as.character(K)]][[sprintf('%s_%s', group1, group2)]]
              
              # Create binary matrices for gene ranks
              cv.train <- train[which(train.label == group1 | train.label == group2),]
              cv.train.label <- factor(train.label[train.label == group1 | train.label == group2])
              cv.train <- cv.train[ ,ktsp$geneX] - cv.train[ ,ktsp$geneY]
              cv.train <- as.data.frame(ifelse(cv.train > 0, "1", "0"))
              colnames(cv.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
              cv.train <- cbind.data.frame(cv.train, type=cv.train.label)
              
              cv.valid <- valid[which(valid.label == group1 | valid.label == group2),]
              cv.valid.label <- factor(valid.label[valid.label == group1 | valid.label == group2])
              cv.valid <- cv.valid[ ,ktsp$geneX] - cv.valid[ ,ktsp$geneY]
              cv.valid <- as.data.frame(ifelse(cv.valid > 0, "1", "0"))
              colnames(cv.valid) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
              cv.valid <- cbind.data.frame(cv.valid, type=cv.valid.label)
              for (p in names(cv.valid)) {
                if (p != "type") {
                  levels(cv.valid[[p]]) <- c("0", "1") -> levels(cv.valid[[p]])
                }
              } 
              
              # Build random forest for pairwise groups
              rf.model <- randomForest(type~., data = cv.train, ntree = 501, replace = F) # tuned later
              #cv.result[[i]][['RFmodels']][[sprintf('%s_%s', group1, group2)]] <- rf.model
              predictionPrs <- as.data.frame(predict(rf.model, newdata=cv.valid, type='prob'))
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
          predictionScore <- opti.RFs(prs.table, grp.names)
          predictionScore$predict.grp <- colnames(predictionScore)[apply(predictionScore, 1, which.max)] -> predictionLabel
        }
      }
      
      # evaluation
      cm <- confusionMatrix(data=predictionLabel, reference=valid.label)
      cv.result[[i]][['model']] <- rf.model
      cv.result[[i]][['trueLabel']] <- valid.label
      cv.result[[i]][['predictionLabel']] <- predictionLabel
      cv.result[[i]][['predictionScore']] <- predictionScore
      cv.result[[i]][['ConfusionMatrix']] <- cm
    }
    
    # acquire overall result
    all.true.label <- c() -> all.predict.label
    all.predict.score <- data.frame()
    for (f in 1:length(cv.result)) {
      all.true.label <- c(all.true.label, as.character(cv.result[[f]][['trueLabel']]))
      all.predict.label <- c(all.predict.label, as.character(cv.result[[f]][['predictionLabel']]))
      all.predict.score <- rbind(all.predict.score, cv.result[[f]][['predictionScore']])
    }
    
    all.cm <- confusionMatrix(data=factor(all.predict.label), reference=factor(all.true.label))
    all.cv[[as.character(t)]][['cm']] <- all.cm
    all.cv[[as.character(t)]][['trueLabel']] <- all.true.label
    all.cv[[as.character(t)]][['ACC']] <- as.numeric(all.cm$overall['Accuracy'])
    all.cv[[as.character(t)]][['Sens']] <- as.numeric(cm$byClass[ ,1]) -> sens
    all.cv[[as.character(t)]][['Spec']] <- as.numeric(cm$byClass[ ,2]) -> spec
    all.cv[[as.character(t)]][['Youden']] <- sens + spec - 1
    names(all.cv[[as.character(t)]][['Youden']]) <- gsub('Class: ', '', rownames(cm$byClass))
    all.cv[[as.character(t)]][['predictionScore']] <- all.predict.score
    
    all.acc <- c(all.acc, all.cv[[as.character(t)]][['ACC']])
    all.sens <- rbind(all.sens, all.cv[[as.character(t)]][['Sens']])
    all.spec <- rbind(all.spec, all.cv[[as.character(t)]][['Spec']])
    all.youden <- rbind(all.youden, all.cv[[as.character(t)]][['Youden']])
  }
  
  all.youden$meanV <- apply(all.youden, 1, mean)
  all.sens$meanV <- apply(all.sens, 1, mean)
  all.spec$meanV <- apply(all.spec, 1, mean)
  all.cv$cm <- all.cm
  all.cv$all.acc <- all.acc
  all.cv$all.sens <- all.sens
  all.cv$all.spec <- all.spec
  all.cv$all.youden <- all.youden
  print(paste0('AVG Youden: ', round(mean(all.youden$meanV), digits = 3)))
  all.cv
}


