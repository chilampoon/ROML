# Repeated 5-fold cross validation
# Inputed data.train: dataframe, sample by gene

suppressPackageStartupMessages({
  library(limma)
  library(caret)
  library(caTools)
  library(randomForest)
  library(pROC)
})
source('~/ktsp/script/function/cal.TSPscore.R')


bi.cv <- function(time=10, data.train, lb.train, feature.num, baseline=TRUE) {
  if (ncol(data.train) != length(lb.train)) stop('Matrix input should be a dataframe which is gene by sample')
  
  # for final benchmarking
  all.cv <- list()
  final.res <- data.frame(feature.num=feature.num, MEANbac=rep(0, length(feature.num)), SDbac=rep(0, length(feature.num)),
                          sens=rep(0, length(feature.num)), spec=rep(0, length(feature.num)), 
                          MEANyd=rep(0, length(feature.num)), SDyd=rep(0, length(feature.num)))
  
  # repeat n times
  for (t in 1:time) {
    print(paste0('Time = ', t))
    # Create 5 folds
    folds <- createFolds(lb.train, k=5)
    one.res <- data.frame(feature.num=feature.num, MEANbac=rep(0, length(feature.num)), sens=rep(0, length(feature.num)),
                          spec=rep(0, length(feature.num)), MEANyd=rep(0, length(feature.num)))
    perf <- list()
    for (num in feature.num) {perf[[as.character(num)]] <- list()}
    
    # loop each fold
    for (i in 1:5) {
      # 1/5 for validation, 4/5 for training
      train <- data.train[ ,-folds[[i]]]
      train.label <- factor(lb.train[-folds[[i]]])
      valid <- data.train[ ,folds[[i]]]
      valid.label <- factor(lb.train[folds[[i]]])
      valid <- cbind(as.data.frame(t(valid)), type = valid.label)
      
      # loop each feature number
      for (num in feature.num) {
        if (baseline == TRUE) {
          #print('running baseline...')
          ## 1. differential expression for baseline feature selection
          type <- train.label
          design <- model.matrix(~ type)
          
          # fitting - note that it's limma-trend here
          vfit <- lmFit(train, design)
          efit <- eBayes(vfit, trend = T)
          
          # get top table
          tt.train <- topTable(efit, number = Inf, coef = colnames(design)[2])
          tt.train <- tt.train[tt.train$adj.P.Val < 0.05, ]
  
          if (num > nrow(tt.train)) {
            stop(paste0('Feature number ', num, ' is too large, set a smaller one.'))
          } else {
            feature <- rownames(tt.train[c(1:num), ])
          }
          
          sub.train <- cbind(as.data.frame(t(train[match(feature, rownames(train)), ])), type = train.label)
          sub.valid <- valid[ ,match(c(feature, 'type'), colnames(valid))]
          
          # train a baseline model
          rf.model <- randomForest(type~., data=sub.train, ntree=501, replace=F)
        } else {
          #print('running ROML...')
          tmp.tsp <- cal.TSPscore(train, train.label)
          ktsp <- tmp.tsp[c(1:num), ]
          
          # Create binary matrices for gene ranks
          bi.train <- train[match(ktsp$geneX, rownames(train)), ] - train[match(ktsp$geneY, rownames(train)), ]
          bi.train <- as.data.frame(ifelse(bi.train > 0, "1", "0"))
          rownames(bi.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
          bi.train <- cbind(as.data.frame(t(bi.train)), type=train.label)
          
          sub.valid <- valid[match(ktsp$geneX, rownames(valid)), ] - valid[match(ktsp$geneY, rownames(valid)), ]
          sub.valid <- as.data.frame(ifelse(sub.valid > 0, "1", "0"))
          rownames(sub.valid) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
          sub.valid <- cbind(as.data.frame(t(sub.valid)), type=valid.label)
          for (p in names(sub.valid)) { # Ensure the levels of predictors are same in train & validate
            if (p != "type") levels(sub.valid[[p]]) <- c("0", "1") -> levels(bi.train[[p]])
          }
          
          # build model using binary features
          rf.model <- randomForest(type~., data = bi.train, ntree = 501, replace = F) 
        }
        # test on sub.valid
        pred.label <- factor(predict(rf.model, newdata=sub.valid, type='response'), levels = levels(valid.label))
        pred.score <- predict(rf.model, newdata=sub.valid, type='prob')
        
        perf[[as.character(num)]][['all.true.label']] <- c(perf[[as.character(num)]][['all.true.label']], as.character(valid.label))
        perf[[as.character(num)]][['all.pred.label']] <- c( perf[[as.character(num)]][['all.pred.label']], as.character(pred.label))
        perf[[as.character(num)]][['all.pred.score']] <- rbind(perf[[as.character(num)]][['all.pred.score']], pred.score)
      }
    }
    
    # overall result of one-time CV for one feature
    for (num in feature.num) {
      all.true.label <- factor(perf[[as.character(num)]][['all.true.label']], levels=levels(lb.train))
      all.pred.label <- factor(perf[[as.character(num)]][['all.pred.label']], levels=levels(lb.train))
      cm <- confusionMatrix(data=all.pred.label, reference=all.true.label)
      
      one.res[one.res$feature.num==num,]$MEANbac <- cm[["byClass"]][["Balanced Accuracy"]]
      one.res[one.res$feature.num==num,]$sens <- cm[["byClass"]][["Sensitivity"]]
      one.res[one.res$feature.num==num,]$spec <- cm[["byClass"]][["Specificity"]]
      youden <- cm[["byClass"]][["Sensitivity"]] + cm[["byClass"]][["Specificity"]] - 1
      one.res[one.res$feature.num==num,]$MEANyd <- youden
    }
    all.cv[[t]] <- one.res
  }

  # save avg. result of n times in final.res
  for (num in feature.num) {
    for (p in c('MEANbac', 'sens', 'spec', 'MEANyd')) {
      final.res[final.res$feature.num == num,][[p]] <- mean(sapply(c(1:t), function(x) all.cv[[x]][feature.num==num,][[p]]))
    }
    final.res[final.res$feature.num == num,]$SDbac <- ifelse(t!=1, sd(sapply(c(1:t), function(x) all.cv[[x]][feature.num==num,]$MEANbac)), 0)
    final.res[final.res$feature.num == num,]$SDyd <- ifelse(t!=1, sd(sapply(c(1:t), function(x) all.cv[[x]][feature.num==num,]$MEANyd)), 0)
  }
  print(final.res)
  list(all.cv=all.cv, res.df=final.res)
}
        
