split_statistic <- function(gene.pair, data.train, label, method = "min"){
  gene1 <- data.train[unlist(gene.pair)[1],]
  gene2 <- data.train[unlist(gene.pair)[2],]
  delta <- as.numeric(gene1 - gene2)
  class <- levels(label)
  delta.0 <- delta[which(label == class[1])]
  delta.1 <- delta[which(label == class[2])]
  if (method == "min") {
    t0 <- mean(delta.0)/sqrt(var(delta.0)/length(delta.0))
    t1 <- mean(delta.1)/sqrt(var(delta.1)/length(delta.1))
    t <- min(abs(t0), abs(t1))
    res <- c(t0, t1, t)
    names(res) <- paste0(c(as.character(class[1]), as.character(class[2]), "min.abs"), '.t')
  } else if (method == "t.test") {
    res <- t.test(delta.0, delta.1, var.equal = F)
  } else {
    stop ("method has to be min or t.test")
  }
  res
}

