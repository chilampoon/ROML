## Functions for plotting

suppressPackageStartupMessages({
  library(RColorBrewer)
  library(tidyr)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(stringr)
})



# Create a color vector containing up to 74 colour hexcodes
getColors <- function(n) {
  col <- brewer.pal.info[brewer.pal.info$category=='qual', ] # get max. 74 colours
  col_vector <- unlist(mapply(brewer.pal, col$maxcolors, rownames(col)))
  ifelse (n > length(col_vector), 
          vec <- sample(col_vector, n, replace=T),
          vec <- sample(col_vector, n, replace=F)
  )
  vec
}


# Density plot
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


# Plot the distribution of two dataframe 
# Input countMat: gene by sample
dataDis <- function(data.list, labels) {
  # Transform all values into one df
  all.long <- data.frame()
  for (i in 1:length(data.list)) {
    data <- as.data.frame(data.list[[i]])
    long.data <- data %>% gather(ID)
    long.data$Data <- labels[i]
    all.long <- rbind(all.long, long.data)
  }
  all.long$Data <- factor(all.long$Data, levels = labels)
  
  # calculate mean for each category
  mu <- all.long %>%
    select(value, Data) %>%
    group_by(Data) %>%
    summarise(mean.value = mean(value, na.rm = T))
  
  ggplot(all.long, aes(x=value, fill=Data)) + 
    #geom_histogram(aes(y=..density..), alpha=.5, binwidth=0.6) +
    geom_density(alpha=0.4) +
    geom_vline(data=mu, aes(xintercept=mean.value, color=Data), linetype='dashed', alpha= 0.5) + 
    scale_fill_brewer(palette="Pastel1") + 
    labs(x='Value', y='Density') +
    theme_bw()
}


add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}


# plot(NA,xlim=c(1, ncol(brca.train)), ylim=c(5,15), xlab="sample", ylab="expression")
# for(i in c('MYH11', 'SPAG5')) {
#   lines(x=1:ncol(brca.train), y=brca.train[i,], col=getColors(2))
# }

