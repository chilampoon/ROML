# Quantile normalization of training and testing data
# The distribution of testing data will follow the one of training data
# Only need to do this when input the exact expression values into ML models

# train.mat, test.mat: using dataframe is better

quantile.norm <- function(train.mat, test.mat) {
  library(preprocessCore, quietly = T)
  
  target <- normalize.quantiles.determine.target(train.mat)
  nm.test.mat <- normalize.quantiles.use.target(as.matrix(test.mat), target)
  
  # Use the original rownames & colnames
  nm.data.test <- as.data.frame(nm.data.test)
  rownames(nm.test.mat) <- rownames(test.mat)
  colnames(nm.test.mat) <- colnames(test.mat)
  
  nm.test.mat
}
