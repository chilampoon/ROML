# Quantile normalization of training and testing data
# The distribution of testing data will follow the one of training data
# Only need to do this when input the exact expression values into ML models

# train.mat, test.mat: numeric matrix only

quantile.norm <- function(train.mat, test.mat) {
  library(preprocessCore, quietly = T)
  
  target <- normalize.quantiles.determine.target(train.mat)
  nm.test.mat <- normalize.quantiles.use.target(as.matrix(test.mat), target)
  
  # Use the original rownames & colnames
  nm.test.mat <- as.data.frame(nm.test.mat)
  rownames(nm.test.mat) <- rownames(test.mat)
  colnames(nm.test.mat) <- colnames(test.mat)
  
  nm.test.mat
}
