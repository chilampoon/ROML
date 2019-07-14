# Function to extract RNA-seq parameters from a real dataset
# Copied from from polyester package: https://github.com/alyssafrazee/polyester/blob/master/R/get_params.R

get.params <- function(counts, threshold=NULL) {
  
  if (!is.null(threshold)) {
    rowm <- rowMeans(counts)
    index1 <- which(rowm > threshold)
    counts <- counts[index1,]
  }
  
  nsamples <- dim(counts)[2]
  counts0 <- counts==0
  nn0 <- rowSums(!counts0)
  if (any(nn0 == 1)) {
    # need more than 1 nonzero count to estimate variance
    counts <- counts[nn0 > 1, ]
    nn0 <- nn0[nn0 > 1]
    counts0 <- counts==0
  }
  mu <- rowSums((!counts0)*counts)/nn0 # The estimated negative binomial mean by method of moments for the non-zero counts
  s2 <- rowSums((!counts0)*(counts - mu)^2)/(nn0-1)
  size <- mu^2/(s2-mu + 0.0001)
  size <- ifelse(size > 0, size, min(size[size > 0])) # The estimated negative binomial size by method of moments for the non-zero counts
  p0 <- (nsamples-nn0)/nsamples # A vector of probabilities that the count will be zero, one for each gene/transcript.
  
  lsize <- log(size)
  lmu <- log(mu + 0.0001)
  fit <- smooth.spline(lsize ~ lmu) # A fit relating log mean to log size for use in simulating new data
  list(p0=p0, mu=mu, size=size, fit=fit)
}
