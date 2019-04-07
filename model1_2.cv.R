tsp.rf.cv <- function(nfold, Ks, labels, dataframe, kpairs) {
  # Create folds
  folds <- createFolds(labels, k = nfold)
  results.list <- list()
  
  # Outer loop for cross-validation
  # Inner loop for different Ks
  for (i in 1:5) {
    results.list[[i]] <- list()
    # 1/5 of the dataset as validation, remaining data for training
    cv.train <- dataframe[-folds[[i]], ]
    cv.train.labels <- labels[-folds[[i]]]
    cv.validate <- dataframe[folds[[i]], ]
    cv.validate.labels <- labels[folds[[i]]]
    
    for (K in Ks) {
      results.list[[i]][[as.character(K)]] <- list()
      all.tsp <- kpairs[[as.character(K)]]
      results.list[[i]][[as.character(K)]][['kpairs']] <- all.tsp
      
      
      # Create binary matrices for gene ranks
      bi.cv.train <- cv.train[, all.tsp$geneX] - cv.train[, all.tsp$geneY]
      bi.cv.train <- as.data.frame(ifelse(bi.cv.train > 0, "1", "0"))
      colnames(bi.cv.train) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
      bi.cv.train <- cbind.data.frame(bi.cv.train, type=cv.train.labels)
      
      bi.cv.validate <- cv.validate[, all.tsp$geneX] - cv.validate[, all.tsp$geneY]
      bi.cv.validate <- as.data.frame(ifelse(bi.cv.validate > 0, "1", "0"))
      colnames(bi.cv.validate) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
      bi.cv.validate <- cbind.data.frame(bi.cv.validate, type=cv.validate.labels)
      for (p in names(bi.cv.validate)) {
        if (p != "type") {
          levels(bi.cv.validate[[p]]) <- c("0", "1") -> levels(bi.cv.train[[p]])
        }
      } # Ensure the levels of predictors are same in train & validate
      
      # Build random forest
      rf.model <- randomForest(type~., data = bi.cv.train, ntree = 500, replace = F) # tuned later
      predictionLabels <- predict(rf.model, newdata=bi.cv.validate, type='response')
      predictionPrs <- predict(rf.model, newdata=bi.cv.validate, type='prob')
      
      # Store performances
      cm <- confusionMatrix(data=as.factor(predictionLabels), reference=as.factor(cv.validate.labels))
      sensitivity <- as.numeric(cm$byClass[, 1])
      specificity <- as.numeric(cm$byClass[, 2])
      youden <- sensitivity + specificity - 1
      names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
      results.list[[i]][[as.character(K)]][['trueLabels']] <- cv.validate.labels
      results.list[[i]][[as.character(K)]][['predictionLabels']] <- predictionLabels
      results.list[[i]][[as.character(K)]][['predictionScores']] <- predictionPrs
      results.list[[i]][[as.character(K)]][['ConfusionMatrix']] <- cm
      results.list[[i]][[as.character(K)]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
      results.list[[i]][[as.character(K)]][['Youden']] <- youden
      results.list[[i]][[as.character(K)]][['model']] <- rf.model
      
      print(paste0("Fold ", i, " & K=", K,  " ACC:", round(results.list[[i]][[as.character(K)]][['ACC']], 5)))
    }
  }
  results.list
}
