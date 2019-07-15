# 5-fold cross validation

# Inputed data.train: dataframe, sample by gene

rf.CV <- function(data.train, lb.train) {
  results.list <- list()
  if (nrow(data.train) != length(lb.train)) stop('Matrix input should be a dataframe which is sample by gene.')
  data.train <- cbind.data.frame(data.train, type=lb.train)
  folds <- createFolds(lb.train, k=5)
  
  for (i in 1:5) {
    results.list[[as.character(i)]] <- list()
    
    # 1/5 for validation, 4/5 for training
    cv.train <- data.train[-folds[[i]], ]
    cv.train.labels <- as.factor(lb.train[-folds[[i]]])
    cv.valid <- data.train[folds[[i]], ]
    cv.valid.labels <- as.factor(lb.train[folds[[i]]])
    
    # Build the model
    rf.model <- randomForest(type~., data=cv.train, ntree=500, replace=F)
    predictionLabels <- predict(rf.model, newdata=cv.valid, type='response')
    predictionScores <- predict(rf.model, newdata=cv.valid, type='prob')
    
    # Save performances
    cm <- confusionMatrix(data=as.factor(predictionLabels), reference=cv.valid.labels)
    sens <- as.numeric(cm$byClass['Sensitivity'])
    spec <- as.numeric(cm$byClass['Specificity'])
    youden <- sens + spec - 1
    
    results.list[[i]][['trueLabels']] <- cv.valid.labels
    results.list[[i]][['predictionLabels']] <- predictionLabels
    results.list[[i]][['predictionScores']] <- predictionScores
    results.list[[i]][['ConfusionMatrix']] <- cm
    results.list[[i]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
    results.list[[i]][['Youden']] <- youden
    results.list[[i]][['AUC']] <- auc(roc(response=as.vector(cv.valid.labels), predictor=as.vector(predictionScores[ ,1]), quiet = T))
    results.list[[i]][['model']] <- rf.model
    print(paste0("Fold ", i, "-ACC:", round(results.list[[i]][['ACC']], 3)))
  }
  results.list
}
