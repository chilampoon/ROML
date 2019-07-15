# Test the ML model with optimized parameters

testing <- function(data.train, data.test, labels.train, labels.test) {
  if (nrow(data.train) != length(lb.train) | nrow(data.test) != length(labels.test)) 
    stop('Matrix input should be a dataframe which is sample by gene.')
  
  # Train a model using all training data (maybe add parameters options later)
  data.train <- cbind.data.frame(data.train, type=as.factor(labels.train))
  full.model <- randomForest(type~., data = data.train, ntree = 500, replace = F)
  
  # Testing on independent data
  test.results <- list()
  data.test <- cbind.data.frame(data.test, type=labels.test)
  
  predictionLabels.test <- predict(model, newdata = data.test, type = 'response')
  predictionLabels.test <- factor(predictionLabels.test, levels = levels(labels.test))
  predictionScores.test <- predict(model, newdata = data.test, type = 'prob')
  
  # Save performances
  cm <- confusionMatrix(data=predictionLabels.test, reference=labels.test)
  sens <- as.numeric(cm$byClass['Sensitivity'])
  spec <- as.numeric(cm$byClass['Specificity'])
  youden <- sens + spec - 1
  #names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
  
  test.results[['model']] <- full.model
  test.results[['ACC']] <- as.numeric(cm$overall['Accuracy'])
  test.results[['predictionLabels']] <- predictionLabels.test
  test.results[['predictionScores']] <- predictionScores.test
  test.results[['cm']] <- cm
  test.results[['AUC']] <- auc(roc(response=as.vector(labels.test), predictor=as.vector(predictionScores.test[,1]), quiet = T))
  test.results[['Youden']] <- youden
  print(paste0("Testing ACC:", round(test.results$ACC, 3)))
  test.results
}
