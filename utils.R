####### Functions
# 1. getColors()
# 2. Plot ACC & AUC performances using box plot
# 3. Get overall AUC or ACC from cv.results
# 4. Plot accuracy
# 5. Plot Youden index
# 6. low means filter
# 7. Plot prediction scores distribution

suppressPackageStartupMessages({
  library(RColorBrewer)
  library(tidyr)
  library(reshape2)
})

source('~/GoogleDrive/pjs/kTSP/code/simu.multi.R')
source('~/GoogleDrive/pjs/kTSP/code/get.ktsp.R')
source('~/GoogleDrive/pjs/kTSP/code/cv.rfTSP.R')
source('~/GoogleDrive/pjs/kTSP/code/testing.R')



# 1. Create a color vector containing up to 74 colour hexcodes
getColors <- function(n) {
  col <- brewer.pal.info[brewer.pal.info$category=='qual', ] # get max. 74 colours
  col_vector <- unlist(mapply(brewer.pal, col$maxcolors, rownames(col)))
  ifelse (n > length(col_vector), 
          vec <- sample(col_vector, n, replace=T),
          vec <- sample(col_vector, n, replace=F)
  )
  vec
}



# 4. Accuracy
plot.acc <- function(overall, Ks, yrange) {
  acc.df <- data.frame(Ks = factor(names(overall), names(overall)), ACC = rep(0, length(overall)))
  for (k in 1:length(Ks)) {
    K <- Ks[k]
    acc.df[k, ]$ACC <- overall[[as.character(K)]][['ACC']]
  }
  
  ggplot(acc.df, aes(x=Ks, y=ACC, group=1)) + geom_line(size=1.1, color="steelblue") + geom_point() +
    theme_bw() +
    scale_y_continuous(breaks = seq(yrange[1], yrange[2], by=0.1), limits = yrange)
  #coord_cartesian(ylim=yrange) 
  #theme(aspect.ratio = 4/3)
}

# 5. Youden index
plot.youden <- function(overall, Ks, yrange) {
  yd.df <- data.frame()
  for (k in 1:length(Ks)) {
    K <- Ks[k]
    yd.df <- rbind.data.frame(yd.df, data.frame(as.list(overall[[as.character(K)]][["Youden"]])))
  }
  yd.df <- cbind.data.frame(Ks = factor(names(overall), names(overall)), yd.df)
  yd.df <- melt(yd.df, id.vars = "Ks")
  
  ggplot(yd.df, aes(x=Ks, y=value, group=variable)) + 
    geom_line(aes(color=variable)) + geom_point(aes(color=variable)) +
    scale_y_continuous(breaks = seq(yrange[1], yrange[2], by=0.1), limits = yrange) +
    theme_bw() + scale_color_brewer(palette="Set1") + labs(y="Youden index", color = "Subtypes")
  
  # ggplot(yd.df, aes(x=variable, y=value, fill=Ks)) + 
  #   geom_bar(stat="identity", position=position_dodge()) +
  #   theme_minimal() + scale_fill_brewer(palette="Accent") + labs(y="Youden index", x = "Types", fill = "Ks")
}


# 6. Low means (of expression values) filter
# cal.mat: X-sample; Y-genes; percent: how many percents of genes you want
lowMeans.filter <- function(cal.mat, percent) {
  means <- data.frame(means = rowMeans(cal.mat), ranks = rank(rowMeans(cal.mat)))
  cuts <- quantile(means$ranks, probs = seq(0, 1, 0.05))
  
  # Top n% genes (rank() is from the smallest to largest)
  cut.percent <- paste0((100 - percent), "%")
  means <- means[which(means$ranks >= cuts[cut.percent]), ]
  sub.genes <- rownames(means)
  sub.genes
}




# 10. Draw the prediction score distribution

prob.jitt <- function(overall, K, label.df, grp.order) {
  #grp.order <- c('Basal', 'Her2', 'LumA', 'LumB', 'Normal')
  prs.table <- as.data.frame(overall[[as.character(K)]][["predictionScores"]])
  if ('predict.grp' %in% colnames(prs.table)) prs.table <- prs.table[, names(prs.table) != 'predict.grp'] 
  prs.table <- prs.table[, grp.order]
  prs.table <- merge(prs.table, label.df, by=0)
  plot.df <- prs.table[, -1]
  plot.df <- melt(plot.df, id.vars = colnames(label.df))
  colnames(plot.df) <- c("true_label", "Pred_Type", "variable")
  
  ggplot(plot.df, aes(x=Pred_Type, y=variable, color=Pred_Type, shape=Pred_Type)) + 
    geom_jitter(position=position_jitter(0.2), cex=.5) + 
    scale_y_continuous(breaks = seq(0, 1, by=0.25), limits = c(-0.005, 1.005)) +
    facet_wrap(~true_label) + theme_bw() + 
    scale_color_brewer(palette="Dark2") + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(x="Predictive Type", y="Probability") 
}



# 11. Draw the density plot for expression data
dens.plot <- function(table, colVec, yrange) {
  d <- plot(density(table[, 1]), col=colVec[1], 
            lwd=2, las=2, ylim=yrange, main="", xlab="") +
    abline(v=0, lty=3) + title(xlab="expr values") +
    for (i in 2:ncol(table)) {
      den <- density(table[, i])
      lines(den$x, den$y, col=colVec[i], lwd=2)
    } 
  d
}
