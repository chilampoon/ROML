object.rf.cv <- function(nfold, Ks, labels, dataframe, kpairs) {
  # Create folds
  folds <- createFolds(labels, k = nfold)
  results.list <- list()
  grp.names <- as.character(unique(labels)) #concordant with the order in TSP score calculation
  
  # Outer loop for cross-validation
  # Intermediate loop for different Ks
  # Inner loop for different classes
  for (i in 1:5) {
    results.list[[i]] <- list()
    # 1/5 of the dataset as validation, remaining data for training
    cv.train <- dataframe[-folds[[i]], ]
    cv.train.labels <- labels[-folds[[i]]]
    cv.validate <- dataframe[folds[[i]], ]
    cv.validate.labels <- labels[folds[[i]]]
    
    for (K in Ks) {
      print(paste0('Loop', i, ' K = ', K, '...'))
      results.list[[i]][[as.character(K)]] <- list()
      results.list[[i]][[as.character(K)]][['RFmodels']] <- list()

      # Get RF models from each group pair
      print('  Building RF models from each group pair')
      prs.table <- data.frame()
      for (g in 1:(length(grp.names)-1)) {
        for (r in (g+1):length(grp.names)) {
          group1 <- grp.names[g]
          group2 <- grp.names[r]
          ktsp <- kpairs[[sprintf('%s_%s', group1, group2)]][[as.character(K)]]
          
          # Create binary matrices for gene ranks
          bi.cv.train <- cv.train[which(cv.train.labels == group1 | cv.train.labels == group2),]
          bi.cv.train.labels <- factor(cv.train.labels[cv.train.labels == group1 | cv.train.labels == group2])
          bi.cv.train <- bi.cv.train[, ktsp$geneX] - bi.cv.train[, ktsp$geneY]
          bi.cv.train <- as.data.frame(ifelse(bi.cv.train > 0, "1", "0"))
          colnames(bi.cv.train) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
          bi.cv.train <- cbind.data.frame(bi.cv.train, type=bi.cv.train.labels)
          
          bi.cv.validate <- cv.validate[, ktsp$geneX] - cv.validate[, ktsp$geneY]
          bi.cv.validate <- as.data.frame(ifelse(bi.cv.validate > 0, "1", "0"))
          colnames(bi.cv.validate) <- sprintf("%s.%s", ktsp$geneX, ktsp$geneY)
          bi.cv.validate <- cbind.data.frame(bi.cv.validate, type=cv.validate.labels)
          for (p in names(bi.cv.validate)) {
            if (p != "type") {
              levels(bi.cv.validate[[p]]) <- c("0", "1") -> levels(bi.cv.train[[p]])
            }
          } # Ensure the levels of predictors are same in train & validate
          
          # Build random forest for pair-wise groups
          rf.model <- randomForest(type~., data = bi.cv.train, ntree = 500, replace = F) # tuned later
          results.list[[i]][[as.character(K)]][['RFmodels']][[sprintf('%s_%s', group1, group2)]] <- rf.model
          predictionPrs <- as.data.frame(predict(rf.model, newdata=bi.cv.validate, type='prob'))
          predictionPrs <- predictionPrs[, c(group1, group2)]
          colnames(predictionPrs) <- paste0(colnames(predictionPrs), paste0(".", sprintf('%s_%s', group1, group2)))
          
          # Combine & merge probability tables
          if (nrow(prs.table) == 0) {
            prs.table <- rbind.data.frame(prs.table, predictionPrs)
          } else {
            prs.table <- cbind.data.frame(prs.table, predictionPrs)
          }
        }
      }
      
      # Optimize final probabilities using objective function
      opti.RFs <- function(prs.table, grp.names) {
        lb <- rep(0, length(grp.names)) # lower bound
        ub <- rep(1, length(grp.names)) # upper bound
        x0 <- setNames(rep(1/length(grp.names), length(grp.names)), grp.names) # vector with starting values for the optimization
        
        # Objective function
        obj_f <- function(x, pw.p.mat) {
          # x0 Normal   LumA   LumB   Her2  Basal 
          #      0.2    0.2    0.2    0.2    0.2 
          fml <- 0
          for (g in 1:(length(x)-1)) {
            for (r in (g+1):length(x)) {
              fml <- fml + (pw.p.mat[g, r] - x[g]/(x[g] + x[r]))^2 + (pw.p.mat[r, g] - x[r]/(x[g] + x[r]))^2
            }
          }
          sum(fml)
        }
        
        # Rewrite the equality as 1-(P1+P2+P3+P4+P5)=0
        eq_c <- function(x, pw.p.mat) 1 - sum(x) 
      
        # Rewrite the inequality as P1*P2*P3*P4*P5 - 1 <= 0
        ineq_c <- function(x, pw.p.mat) prod(x) - 1 
        
        # Function to transform a row in prs.table to a matrix
        rowToMat <- function(row, grp.names) {
          row.mat <- setNames(data.frame(matrix(ncol = 5, nrow = 5), stringsAsFactors=F), grp.names)
          rownames(row.mat) <- colnames(row.mat)
          for (g in 1:(length(grp.names)-1)) {
            grp1 <- grp.names[g]
            for (r in (g+1):length(grp.names)) {
              grp2 <- grp.names[r]
              grp <- sprintf('%s_%s', grp1, grp2)
              row.mat[g, r] <- row[, sprintf('%s.%s', grp1, grp)]
              row.mat[r, g] <- row[, sprintf('%s.%s', grp2, grp)]
            }
          }
          row.mat
        }
        
        
        # Summarize final probabilities
        final.prs <- setNames(data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors=F), grp.names)
        
        for (n in 1:nrow(prs.table)) {
          pw.p.mat <- rowToMat(prs.table[n, ], grp.names)
          opt <- nloptr(x0=x0, eval_f=obj_f, lb=lb, ub=ub, eval_g_ineq=ineq_c, eval_g_eq=eq_c, 
                        opts=list(algorithm="NLOPT_GN_ISRES", "maxeval" = 150000), pw.p.mat=pw.p.mat)
          final.prs[nrow(final.prs)+1, ] <- opt$solution
        }
        rownames(final.prs) <- rownames(prs.table)
        final.prs
      }
      
      print('  Calculating final probabilities using our obj function')
      final.prs <- opti.RFs(prs.table, grp.names)
      final.prs$predict.grp <- colnames(final.prs)[apply(final.prs, 1, which.max)] -> predictionLabels
      
      # Save performances
      cm <- confusionMatrix(data=as.factor(predictionLabels), reference=as.factor(cv.validate.labels))
      sensitivity <- as.numeric(cm$byClass[, 1])
      specificity <- as.numeric(cm$byClass[, 2])
      youden <- sensitivity + specificity - 1
      names(youden) <- gsub('Class: ', '', rownames(cm$byClass))
      results.list[[i]][[as.character(K)]][['trueLabels']] <- cv.validate.labels
      results.list[[i]][[as.character(K)]][['predictionLabels']] <- predictionLabels
      results.list[[i]][[as.character(K)]][['predictionPrsTable']] <- final.prs
      results.list[[i]][[as.character(K)]][['ConfusionMatrix']] <- cm
      results.list[[i]][[as.character(K)]][['ACC']] <- as.numeric(cm$overall['Accuracy'])
      results.list[[i]][[as.character(K)]][['Youden']] <- youden
    }
  }
  results.list
}
