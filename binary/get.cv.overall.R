# Extract and calculate overall performance of cross-validation
# Input the result list from cross validation funciton like rf.CV

get.cv.overall <- function(cv.results.list) {
  # Prepare objects
  overall.results <- list()
  all.true.labels <- c()
  all.predict.labels <- c()
  all.predict.scores <- c()
  
  for (f in 1:length(cv.results.list)) {
    all.true.labels <- c(all.true.labels, as.character(cv.results.list[[f]][['trueLabels']]))
    all.predict.labels <- c(all.predict.labels, as.character(cv.results.list[[f]][['predictionLabels']]))
    all.predict.scores <- rbind(all.predict.scores, cv.results.list[[f]][['predictionScores']])
  }
  
  all.cm <- confusionMatrix(data=as.factor(all.predict.labels), reference=as.factor(all.true.labels))
  overall.results[['trueLabels']] <- all.true.labels
  overall.results[['ACC']] <- as.numeric(all.cm$overall['Accuracy'])
  overall.results[['Sensitivity']] <- as.numeric(all.cm$byClass['Sensitivity']) -> all.sens
  overall.results[['Specificity']] <- as.numeric(all.cm$byClass['Specificity']) -> all.spec
  all.youden <- all.sens + all.spec - 1
  #names(all.youden) <- gsub('Class: ', '', rownames(all.cm$byClass))
  overall.results[['Youden']] <- all.youden
  overall.results[['predictionScores']] <- all.predict.scores
  overall.results[['cm']] <- all.cm
  overall.results[['AUC']] <- auc(roc(response=as.vector(all.true.labels), predictor=as.vector(all.predict.scores[,1]), quiet = T))
  overall.results
}
