# Test the ML model with optimized parameters

# Train a model using all features
train.bi.model <- function(ori.train, lb.train, featureSet) {
  suppressPackageStartupMessages({
    library(randomForest)
  })
  ori.train <- as.data.frame(ori.train)
  train <- ori.train[match(featureSet, rownames(ori.train)),]
  train <- cbind.data.frame(as.data.frame(t(train)), type=lb.train)
  full.model <- randomForest(type~., data = train, ntree = 501, replace = F)
  full.model
}


# test on independent data
bi.test <- function(data.train, data.test, label.train, label.test, model) {
  suppressPackageStartupMessages({
    library(caret)
    library(caTools)
    library(pROC)
    library(randomForest)
  })
  
  if (nrow(data.train) != length(label.train) | nrow(data.test) != length(label.test)) 
    stop('Matrix input should be a dataframe which is sample by gene.')
  
  
  # Testing on independent data
  test.results <- list()
  data.test <- cbind.data.frame(data.test, type=label.test)
  
  predictionlabel.test <- predict(model, newdata = data.test, type = 'response')
  predictionlabel.test <- factor(predictionlabel.test, levels = levels(label.test))
  predictionScore.test <- predict(model, newdata = data.test, type = 'prob')
  
  # Save performances
  cm <- confusionMatrix(data=predictionlabel.test, reference=label.test)
  sens <- as.numeric(cm$byClass['Sensitivity'])
  spec <- as.numeric(cm$byClass['Specificity'])
  youden <- sens + spec - 1
  
  test.results[['model']] <- model
  test.results[['ACC']] <- as.numeric(cm$overall['Accuracy'])
  test.results[['predictionLabels']] <- predictionlabel.test
  test.results[['predictionScores']] <- predictionScore.test
  test.results[['cm']] <- cm
  test.results[['Youden']] <- youden
  test.results[['sens']] <- sens
  test.results[['spec']] <- spec
  test.results[['AUC']] <- auc(roc(response=as.vector(label.test), predictor=as.vector(predictionScore.test[,1]), quiet = T))
  print(paste0("Testing ACC:", round(test.results$ACC, 3)))
  test.results
}
