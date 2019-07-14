## Low mean gene filter

# cal.mat: matrix or dataframe - gene by sample
# percent: percentage of data you want to keep
# return: retained gene ids

lowMeans.filter <- function(cal.mat, percent) {
  means <- data.frame(means = rowMeans(cal.mat), ranks = rank(rowMeans(cal.mat)))
  cuts <- quantile(means$ranks, probs = seq(0, 1, 0.05))
  
  # Top n% genes (rank() is from the smallest to largest)
  cut.percent <- paste0((100 - percent), "%")
  means <- means[which(means$ranks >= cuts[cut.percent]), ]
  sub.genes <- rownames(means)
  sub.genes
}
