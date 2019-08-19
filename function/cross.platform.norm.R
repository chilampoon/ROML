# Cross platform normalization for testing data
# Target - training data; To be normalized - testing data
# Only need to do this when input the exact expression values into ML models


cross.platform.norm <- function(data.train, data.test, method = c('QN', 'TDM', 'NPN')) {
  suppressPackageStartupMessages({
    library(preprocessCore)
    library(TDM)
    library(data.table)
    library(binr)
    library(scales)
  }) 
  source('~/ktsp/script/function/fsqn.R')
  
  if (method == 'QN') {
    # 1. quantile normalization (input: matrix)
    target <- normalize.quantiles.determine.target(as.matrix(data.train))
    qn.data.test <- normalize.quantiles.use.target(as.matrix(data.test), target)
    
    qn.data.test <- as.data.frame(qn.data.test)
    rownames(qn.data.test) <- rownames(data.test)
    colnames(qn.data.test) <- colnames(data.test)
    return(qn.data.test)
    
  } else if (method == 'TDM') {
    # 2. training distribution matching (input: data.table)
    tdm.ref <- data.table(cbind(rownames(data.train) ,data.train))
    tdm.targ <- data.table(cbind(rownames(data.test) ,data.test))
    colnames(tdm.ref)[1] <- 'gene' -> colnames(tdm.targ)[1]
    tdm.data.test <- tdm_transform(target_data=tdm.targ,
                                   ref_data = tdm.ref,
                                   inv_reference = F,
                                   log_target = F)
    
    tdm.data.test <- tdm.data.test[,-1]
    tdm.data.test <- apply(tdm.data.test, c(1,2), as.numeric)
    tdm.data.test <- as.data.frame(tdm.data.test)
    rownames(tdm.data.test) <- rownames(data.test)
    return(tdm.data.test)
   
  } else if (method == 'FSQN') {
    # 3. feature specific quantile normalization (input: matrix)
    fsqn.test <- t(data.test)
    fsqn.target <- t(data.train)
    fsqn.data.test <- quantileNormalizeByFeature(matrix_to_normalize = fsqn.test,
                                                 target_distribution_matrix = fsqn.target)
    
    fsqn.data.test <- as.data.frame(t(fsqn.data.test))
    rownames(fsqn.data.test) <- rownames(data.test)
    colnames(fsqn.data.test) <- colnames(data.test)
    return(fsqn.data.test)
  }
}
