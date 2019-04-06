object.testing <- function(dataframe, labels, grp.names, opti.Ks, kpairs, final.RF) {
  test.results <- list()
  
  for (k in opti.Ks) {
    test.results[[as.character(k)]] <- list()
    
    prs.table.test <- data.frame()
    for (g in 1:(length(grp.names)-1)) {
      for (r in (g+1):length(grp.names)) {
        group1 <- grp.names[g]
        group2 <- grp.names[r]
        ktsp <- kpairs[[sprintf('%s_%s', group1, group2)]][[as.character(k)]]
    
        # Subset data & binarize
        bi.test <- dataframe[, ktsp$geneX] - dataframe[, ktsp$geneY]
        bi.test <- as.data.frame(ifelse(bi.test > 0, "1", "0"))
        colnames(bi.test) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
        bi.test <- cbind.data.frame(bi.test, type=labels)
        for (p in names(bi.test)) {
          if (p != "type") {
            levels(bi.test[[p]]) <- c("0", "1")
          }
        }
        
        model <- final.RF[[as.character(k)]][[sprintf('%s_%s', group1, group2)]]
        predictionPrs.test <- predict(model, newdata=bi.test, type='prob')
        predictionPrs.test <- predictionPrs.test[, c(group1, group2)]
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
    print('  Calculating final probabilities using our obj function')
    final.prs.test <- opti.RFs(prs.table.test, grp.names)
    final.prs.test$predict.grp <- colnames(final.prs.test)[apply(final.prs.test, 1, which.max)] -> predictionLabels.test
    
    
    # Performances
    cm <- confusionMatrix(data=as.factor(predictionLabels.test), reference=as.factor(labels))
    test.acc <- as.numeric(cm$overall['Accuracy'])
    sensitivity <- as.numeric(cm$byClass[, 1])
    specificity <- as.numeric(cm$byClass[, 2])
    youden <- sensitivity + specificity - 1
    names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
    test.results[[as.character(k)]]$ACC <-  test.acc
    test.results[[as.character(k)]]$Youden <- youden
    test.results[[as.character(k)]]$predictionPrsTable <- prs.table.test
    test.results[[as.character(k)]]$cm <- cm
  }
  print(paste0(k, " - the best ACC:", round(test.results[[as.character(k)]]$ACC, 5)))
  test.results
}
