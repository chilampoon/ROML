# Repeated 5-fold cross validation
# Inputed data.train: dataframe, sample by gene

rf.CV <- function(data.train, lb.train, time=10, K, ...) {
  library(caret, quietly = T)
  library(caTools, quietly = T)
  
  if (nrow(data.train) != length(lb.train)) stop('Matrix input should be a dataframe which is sample by gene.')
  data.train <- cbind.data.frame(data.train, type=lb.train)
  
  # to final benchmarking
  all.acc <- c() -> all.spec -> all.sens -> all.youden
  all.cv <- list()
  
  # repeat 10 times
  for (t in 1:time) {
    all.cv[[as.character(t)]] <- list()
    # Create 5 folds
    folds <- createFolds(lb.train, k=5)
    
    # loop each fold
    cv.result <- list()
    for (i in 1:5) {
      cv.result[[i]] <- list()
      # 1/5 for validation, 4/5 for training
      train <- data.train[-folds[[i]], ]
      train.label <- factor(lb.train[-folds[[i]]])
      valid <- data.train[folds[[i]], ]
      valid.label <- factor(lb.train[folds[[i]]])
      
      #### build the model (in two ways) ####
      if (missing(K)) { ## baseline method
        rf.model <- randomForest(type~., data=train, ntree=500, replace=F)
      } else { ## ROML method
        tmp <- list(...)
        all.tsp <- tmp[[1]][[as.character(K)]]
        
        # Create binary matrices for gene ranks
        bi.train <- train[ ,match(all.tsp$geneX, colnames(train))] - cv.train[ ,match(all.tsp$geneY, colnames(train))]
        bi.train <- as.data.frame(ifelse(bi.train > 0, "1", "0"))
        colnames(bi.train) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
        bi.train <- cbind.data.frame(bi.train, type=train.label)
        
        bi.valid <- valid[ ,match(all.tsp$geneX,colnames(valid))] - cv.validate[ ,match(all.tsp$geneY,colnames(valid))]
        bi.valid <- as.data.frame(ifelse(bi.valid > 0, "1", "0"))
        colnames(bi.valid) <- sprintf("%s.%s", all.tsp$geneX, all.tsp$geneY)
        bi.valid <- cbind.data.frame(bi.valid, type=valid.label)
        for (p in names(bi.valid)) { # Ensure the levels of predictors are same in train & validate
          if (p != "type") levels(bi.valid[[p]]) <- c("0", "1") -> levels(bi.train[[p]])
        }
        # build model using binary features
        rf.model <- randomForest(type~., data = bi.train, ntree = 500, replace = F) 
        valid <- bi.valid
      }
      
      predictionLabel <- predict(rf.model, newdata=valid, type='response')
      predictionLabel <- factor(predictionLabel, levels = levels(train.label))
      predictionScore <- predict(rf.model, newdata=valid, type='prob')
      
      # evaluation
      cm <- confusionMatrix(data=predictionLabel, reference=valid.label)
      cv.result[[i]][['model']] <- rf.model
      cv.result[[i]][['trueLabel']] <- valid.label
      cv.result[[i]][['predictionLabel']] <- predictionLabel
      cv.result[[i]][['predictionScore']] <- predictionScore
      cv.result[[i]][['ConfusionMatrix']] <- cm
      #print(paste0("Fold ", i, "-ACC:", round(cv.result[[i]][['ACC']], 3)))
    }
    
    # acquire overall result
    all.true.label <- all.predict.label <- c() -> all.predict.score
    for (f in 1:length(cv.result)) {
      all.true.label <- c(all.true.label, as.character(cv.result[[f]][['trueLabel']]))
      all.predict.label <- c(all.predict.label, as.character(cv.result[[f]][['predictionLabel']]))
      all.predict.score <- rbind(all.predict.score, cv.result[[f]][['predictionScore']])
    }
    
    all.cm <- confusionMatrix(data=factor(all.predict.label), reference=factor(all.true.label))
    all.cv[[as.character(t)]][['cm']] <- all.cm
    all.cv[[as.character(t)]][['trueLabel']] <- all.true.label
    all.cv[[as.character(t)]][['ACC']] <- as.numeric(all.cm$overall['Accuracy'])
    all.cv[[as.character(t)]][['Sens']] <- as.numeric(all.cm$byClass['Sensitivity']) -> all.sens
    all.cv[[as.character(t)]][['Spec']] <- as.numeric(all.cm$byClass['Specificity']) -> all.spec
    all.cv[[as.character(t)]][['Youden']] <- all.sens + all.spec - 1
    all.cv[[as.character(t)]][['predictionScore']] <- all.predict.score
    all.cv[[as.character(t)]][['AUC']] <- auc(roc(response=as.vector(all.true.label), predictor=as.vector(all.predict.score[,1]), quiet = T))
    
    all.acc <- c(all.acc, all.cv[[as.character(t)]][['ACC']])
    all.sens <- c(all.sens, all.cv[[as.character(t)]][['Sens']])
    all.spec <- c(all.spec, all.cv[[as.character(t)]][['Spec']])
    all.youden <- c(all.youden, all.cv[[as.character(t)]][['Youden']])
  }
  all.cv$avgACC <- mean(all.acc)
  all.cv$avgSens <- mean(all.sens)
  all.cv$avgSpec <- mean(all.spec)
  all.cv$avgYouden <- mean(all.youden)
  print(paste0('AVG ACC: ', all.cv$avgACC))
  all.cv
}
