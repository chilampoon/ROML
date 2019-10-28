source('~/ktsp/script/function/objective.fun.R')

suppressPackageStartupMessages({
  library(randomForest)
  library(caret)
  library(caTools)
  library(pROC)
})
  

# Train a model using whole training set
multi.model <- function(dataframe, label, K, grp.names, kpairs, no.model) {
  if (nrow(dataframe) != length(label)) stop('Matrix input should be a dataframe which is sample by gene.')
  
  if (missing(K)) {
    # baseline
    dataframe <- cbind.data.frame(dataframe, type=label)
    model <- randomForest(type~., data = dataframe, ntree = 501, replace = F)
  } else {
    if (no.model != 3) { 
      # For model 1 & 2
      ktsp <- kpairs[[as.character(K)]]
      bi.train <- dataframe[ ,ktsp$geneX] - dataframe[ ,ktsp$geneY]
      bi.train <- as.data.frame(ifelse(bi.train > 0, "1", "0"))
      colnames(bi.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
      bi.train <- cbind.data.frame(bi.train, type=label)
      for (p in names(bi.train)) {
        if (p != "type") levels(bi.train[[p]]) <- c("0", "1") 
      }
      model <- randomForest(type~., data = bi.train, ntree = 501, replace = F) 
    } else {
      # For model 3
      model <- list()
      for (g in 1:(length(grp.names)-1)) {
        for (r in (g+1):length(grp.names)) {
          group1 <- grp.names[g]
          group2 <- grp.names[r]
          ktsp <- kpairs[[as.character(K)]][[sprintf('%s_%s', group1, group2)]]
          bi.train <- dataframe[which(label == group1 | label == group2),]
          bi.label <- factor(label[label == group1 | label == group2])
          bi.train <- bi.train[ ,ktsp$geneX] - bi.train[ ,ktsp$geneY]
          bi.train <- as.data.frame(ifelse(bi.train > 0, "1", "0"))
          colnames(bi.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
          bi.train <- cbind.data.frame(bi.train, type=bi.label)
          for (p in names(bi.train)) {
            if (p != "type") levels(bi.train[[p]]) <- c("0", "1")
          }
          
          # Build random forest for pair-wise groups
          rf.model <- randomForest(type~., data = bi.train, ntree = 501, replace = F) # tuned later
          model[[sprintf('%s_%s', group1, group2)]] <- rf.model
        }
      }
    }
  }
  model
}


# testing
multi.rf.test <- function(testdata, testlabel, model, K, grp.names, kpairs, no.model) {
  if (nrow(testdata) != length(testlabel)) stop('Matrix input should be a dataframe which is sample by gene.')
  test.res <- list()
  
  if (missing(K)) {
    # baseline method
    testdata <- cbind.data.frame(testdata, type=testlabel)
    
    predictionlabel.test <- predict(model, newdata = testdata, type = 'response')
    predictionScore.test <- predict(model, newdata = testdata, type = 'prob')
  } else {
    # ROML methods
    if (no.model != 3) { 
      # For model 1 & 2
      all.tsp <- kpairs[[as.character(K)]]
      bi.test <- testdata[ ,all.tsp$geneX] - testdata[ ,all.tsp$geneY]
      bi.test <- as.data.frame(ifelse(bi.test > 0, "1", "0"))
      colnames(bi.test) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
      bi.test <- cbind.data.frame(bi.test, type=testlabel)
      for (p in names(bi.test)) {
        if (p != "type") levels(bi.test[[p]]) <- c("0", "1") 
      }
      
      predictionlabel.test <- predict(model, newdata = bi.test, type = 'response')
      predictionScore.test <- predict(model, newdata = bi.test, type = 'prob')
    } else { # For model 3
      prs.table.test <- data.frame()
      for (g in 1:(length(grp.names)-1)) {
        for (r in (g+1):length(grp.names)) {
          group1 <- grp.names[g]
          group2 <- grp.names[r]
          ktsp <- kpairs[[sprintf('%s_%s', group1, group2)]]
          
          # Subset data & binarize
          bi.test <- testdata[which(testlabel == group1 | testlabel == group2),]
          new.label <- factor(testlabel[testlabel == group1 | testlabel == group2])
          bi.test <- bi.test[ ,ktsp$geneX] - bi.test[ ,ktsp$geneY]
          bi.test <- as.data.frame(ifelse(bi.test > 0, "1", "0"))
          colnames(bi.test) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
          bi.test <- cbind.data.frame(bi.test, type=new.label)
          for (p in names(bi.test)) if (p != "type") levels(bi.test[[p]]) <- c("0", "1")
          
          model3 <- model[[sprintf('%s_%s', group1, group2)]]
          predictionPrs.test <- predict(model3, newdata=bi.test, type='prob')
          colnames(predictionPrs.test) <- paste0(colnames(predictionPrs.test), paste0(".", sprintf('%s_%s', group1, group2)))
          
          # Combine & merge probability tables
          if (nrow(prs.table.test) == 0) {
            prs.table.test <- rbind.data.frame(prs.table.test, predictionPrs.test)
          } else {
            prs.table.test <- cbind.data.frame(prs.table.test, predictionPrs.test)
          }
        }
      }
      
      # Calculate final probabilities using objective function
      predictionScore.test <- opti.RFs(prs.table.test, grp.names)
      predictionScore.test$predict.grp <- colnames(predictionScore.test)[apply(predictionScore.test, 1, which.max)] -> predictionlabel.test
    }
  }
    
  # Performances
  predictionlabel.test <- factor(predictionlabel.test, levels = levels(testlabel))
  cm <- confusionMatrix(data=predictionlabel.test, reference=testlabel)
  sensitivity <- as.numeric(cm$byClass[,1])
  specificity <- as.numeric(cm$byClass[,2])
  youden <- sensitivity + specificity - 1
  names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
  test.res[['ACC']] <- as.numeric(cm$overall['Accuracy'])
  test.res[['avg.sens']] <- mean(sensitivity)
  test.res[['avg.spec']] <- mean(specificity)
  test.res[['predictionlabel']] <- predictionlabel.test
  test.res[['predictionScore']] <- predictionScore.test
  test.res[['cm']] <- cm
  test.res[['Youden']] <- youden
  test.res[['avg.yd']] <- mean(youden)
  print(paste0("AVG YOUDEN:", round(test.res$avg.yd, 3)))
  test.res
}

