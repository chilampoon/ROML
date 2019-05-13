# Train a model using whole training set
trainModel <- function(dataframe, labels, Ks, grp.names, kpairs, no.model, baseline=FALSE) {
  if (baseline == TRUE) {
    dataframe <- cbind.data.frame(dataframe, type=labels)
    model <- randomForest(type~., data = dataframe, ntree = 500, replace = F)
  } else {
    model <- list()
    for (K in Ks) {
      model[[as.character(K)]] <- list()
      
      if (no.model != 3) { # For model 1 & 2
        ktsp <- kpairs[[as.character(no.model)]][[as.character(K)]]
        bi.train <- dataframe[, ktsp$geneX] - dataframe[, ktsp$geneY]
        bi.train <- as.data.frame(ifelse(bi.train > 0, "1", "0"))
        colnames(bi.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
        bi.train <- cbind.data.frame(bi.train, type=labels)
        for (p in names(bi.train)) {
          if (p != "type") levels(bi.train[[p]]) <- c("0", "1") 
        }
        model[[as.character(K)]] <- randomForest(type~., data = bi.train, ntree = 500, replace = F) 
      } else { # For model 3
        model[[as.character(K)]] <- list()
        for (g in 1:(length(grp.names)-1)) {
          for (r in (g+1):length(grp.names)) {
            group1 <- grp.names[g]
            group2 <- grp.names[r]
            ktsp <- kpairs[[as.character(no.model)]][[as.character(K)]][[sprintf('%s_%s', group1, group2)]]
            bi.train <- dataframe[which(labels == group1 | labels == group2),]
            bi.labels <- factor(labels[labels == group1 | labels == group2])
            bi.train <- bi.train[ ,ktsp$geneX] - bi.train[ ,ktsp$geneY]
            bi.train <- as.data.frame(ifelse(bi.train > 0, "1", "0"))
            colnames(bi.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
            bi.train <- cbind.data.frame(bi.train, type=bi.labels)
            for (p in names(bi.train)) {
              if (p != "type") levels(bi.train[[p]]) <- c("0", "1")
            }
            
            # Build random forest for pair-wise groups
            rf.model <- randomForest(type~., data = bi.train, ntree = 500, replace = F) # tuned later
            model[[as.character(K)]][[sprintf('%s_%s', group1, group2)]] <- rf.model
          }
        }
      }
    }
  }
  model
}

testing <- function(testdata, testLabels, Ks, grp.names, kpairs, full.model, no.model, baseline=FALSE) {
  test.results <- list()
  if (baseline == TRUE) {
    testdata <- cbind.data.frame(testdata, type=testLabels)
    
    predictionLabels.test <- predict(full.model, newdata = testdata, type = 'response')
    predictionLabels.test <- factor(predictionLabels.test, levels = levels(testLabels))
    predictionScores.test <- predict(full.model, newdata = testdata, type = 'prob')
    
    cm <- confusionMatrix(data=predictionLabels.test, reference=testLabels)
    sensitivity <- as.numeric(cm$byClass[, 1])
    specificity <- as.numeric(cm$byClass[, 2])
    youden <- sensitivity + specificity - 1
    names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
    test.acc <- as.numeric(cm$overall['Accuracy'])
    
    test.results[['ACC']] <- test.acc
    test.results[['avg.sens']] <- sum(sensitivity)/length(sensitivity)
    test.results[['avg.spec']] <- sum(specificity)/length(specificity)
    test.results[['predictionLabels']] <- predictionLabels.test
    test.results[['predictionScores']] <- predictionScores.test
    test.results[['cm']] <- cm
    test.results[['Youden']] <- youden
    test.results[['avg.yd']] <- sum(youden)/length(youden)
    print(paste0("ACC:", round(test.results$ACC, 5)))
  } else {
    for (k in Ks) {
      test.results[[as.character(k)]] <- list()
      
      if (no.model != 3) { # For model 1 & 2
        # Subset data & binarize
        all.tsp <- kpairs[[as.character(no.model)]][[as.character(k)]]
        bi.test <- testdata[ ,all.tsp$geneX] - testdata[ ,all.tsp$geneY]
        bi.test <- as.data.frame(ifelse(bi.test > 0, "1", "0"))
        colnames(bi.test) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
        bi.test <- cbind.data.frame(bi.test, type=testLabels)
        for (p in names(bi.test)) {
          if (p != "type") levels(bi.test[[p]]) <- c("0", "1") 
        }
        
        model <- full.model[[as.character(k)]]
        predictionLabels.test <- predict(model, newdata = bi.test, type='response')
        predictionPrs.test <- predict(model, newdata = bi.test, type='prob')
      } else { # For model 3
        prs.table.test <- data.frame()
        for (g in 1:(length(grp.names)-1)) {
          for (r in (g+1):length(grp.names)) {
            group1 <- grp.names[g]
            group2 <- grp.names[r]
            ktsp <- kpairs[[as.character(no.model)]][[as.character(k)]][[sprintf('%s_%s', group1, group2)]]
            
            # Subset data & binarize
            bi.test <- testdata[, ktsp$geneX] - testdata[, ktsp$geneY]
            bi.test <- as.data.frame(ifelse(bi.test > 0, "1", "0"))
            colnames(bi.test) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
            bi.test <- cbind.data.frame(bi.test, type=testLabels)
            for (p in names(bi.test)) {
              if (p != "type") levels(bi.test[[p]]) <- c("0", "1")
            }
            
            model <- full.model[[as.character(k)]][[sprintf('%s_%s', group1, group2)]]
            predictionPrs.test <- predict(model, newdata=bi.test, type='prob')
            predictionPrs.test <- predictionPrs.test[ ,c(group1, group2)]
            colnames(predictionPrs.test) <- paste0(colnames(predictionPrs.test), paste0(".", sprintf('%s_%s', group1, group2)))
            
            # Combine & merge probability tables
            if (nrow(prs.table.test) == 0) {
              prs.table.test <- rbind.data.frame(prs.table.test, predictionPrs.test)
            } else {
              prs.table.test <- cbind.data.frame(prs.table.test, predictionPrs.test)
            }
          }
        }
        
        # Optimize final probabilities using objective function
        #print('  Calculating final probs using our obj. func.')
        predictionPrs.test <- opti.RFs(prs.table.test, grp.names)
        predictionPrs.test$predict.grp <- colnames(predictionPrs.test)[apply(predictionPrs.test, 1, which.max)] -> predictionLabels.test
      }
      
      # Performances
      predictionLabels.test <- factor(predictionLabels.test, levels = levels(testLabels))
      cm <- confusionMatrix(data=predictionLabels.test, reference=testLabels)
      sensitivity <- as.numeric(cm$byClass[, 1])
      specificity <- as.numeric(cm$byClass[, 2])
      youden <- sensitivity + specificity - 1
      names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
      test.results[[as.character(k)]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
      test.results[[as.character(k)]][['avg.sens']] <- sum(sensitivity)/length(sensitivity)
      test.results[[as.character(k)]][['avg.spec']] <- sum(specificity)/length(specificity)
      test.results[[as.character(k)]][['predictionLabels']] <- predictionLabels.test
      test.results[[as.character(k)]][['predictionScores']] <- predictionPrs.test
      test.results[[as.character(k)]][['cm']] <- cm
      test.results[[as.character(k)]][['Youden']] <- youden
      test.results[[as.character(k)]][['avg.yd']] <- sum(youden)/length(youden)
      print(paste0(k, " - ACC:", round(test.results[[as.character(k)]]$ACC, 5)))
    }
  }
  test.results
}