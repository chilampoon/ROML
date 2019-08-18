# Test the ML model with optimized parameters
suppressPackageStartupMessages({
  library(caret)
  library(caTools)
  library(pROC)
  library(randomForest)
})

# Train a model using selected features
train.bi.model <- function(ori.train, lb.train, featureSet, baseline=TRUE) {
  if (ncol(ori.train) != length(lb.train)) stop('Matrix input should be a gene by sample dataframe')
  
  ori.train <- as.data.frame(ori.train)
  if (baseline == TRUE) {
    # baseline model
    train <- ori.train[match(featureSet, rownames(ori.train)),]
    train <- cbind.data.frame(as.data.frame(t(train)), type=lb.train)
  } else {
    # ROML model (featureSet = kpairs dataframe here)
    train <- ori.train[match(featureSet$geneX, rownames(ori.train)), ] - ori.train[match(featureSet$geneY, rownames(ori.train)), ]
    train <- as.data.frame(ifelse(train > 0, "1", "0"))
    rownames(train) <- sprintf("%s.%s", featureSet$geneX, featureSet$geneY)
    train <- cbind(as.data.frame(t(train)), type=lb.train)
    for (p in names(train)) { 
      if (p != "type") levels(train[[p]]) <- c("0", "1") 
    }
  }
  full.model <- randomForest(type~., data = train, ntree = 501, replace = F)
  full.model
}


# test on independent data
bi.test <- function(data.test, label.test, model, kpairs) {
  if (ncol(data.test) != length(label.test)) stop('Matrix input should be a gene by sample dataframe')
  
  test.results <- list()
  if (missing(kpairs)) {
    # baseline method
    final.test <- cbind(as.data.frame(t(data.test)), type=label.test)
  } else {
    # ROML method
    final.test <- data.test[match(kpairs$geneX, rownames(data.test)), ] - data.test[match(kpairs$geneY, rownames(data.test)), ]
    final.test <- as.data.frame(ifelse(final.test > 0, "1", "0"))
    rownames(final.test) <- sprintf("%s.%s", kpairs$geneX, kpairs$geneY)
    final.test <- cbind(as.data.frame(t(final.test)), type=label.test)
    for (p in names(final.test)) {
      if (p != "type") levels(final.test[[p]]) <- c("0", "1")
    }
  }
  
  predictionLabel <- factor(predict(model, newdata = final.test, type = 'response'), levels = levels(label.test))
  predictionScore <- predict(model, newdata = final.test, type = 'prob')
  
  # Save performances
  cm <- confusionMatrix(data=predictionLabel, reference=label.test)
  sens <- as.numeric(cm$byClass['Sensitivity'])
  spec <- as.numeric(cm$byClass['Specificity'])
  youden <- sens + spec - 1
  
  test.results[['model']] <- model
  test.results[['ACC']] <- as.numeric(cm$overall['Accuracy'])
  test.results[['predictionLabel']] <- predictionLabel
  test.results[['predictionScore']] <- predictionScore
  test.results[['cm']] <- cm
  test.results[['Youden']] <- youden
  test.results[['sens']] <- sens
  test.results[['spec']] <- spec
  test.results[['AUC']] <- auc(roc(response=as.vector(label.test), predictor=as.vector(predictionScore[,1]), quiet = T))
  print(paste0("Testing ACC:", round(test.results$ACC, 3)))
  test.results
}
