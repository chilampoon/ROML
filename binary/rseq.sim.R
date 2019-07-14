# Function to generate psudo RNA-seq read counts
# Modified from funciton sim.counts https://github.com/cran/ssizeRNA/blob/master/R/sim.counts.R

rseq.sim <- function(nGenes = 10000, pi0 = 0.8, m, mu, disp, fc, 
                       up = 0.5, replace = TRUE) {
  arg <- list(nGenes = nGenes,
              pi0 = pi0,
              group = rep(c(1, 2), each = m))
  
  ## expected DEG amount
  nDE <- round(nGenes * pi0)
  DE <- nGenes - nDE
  
  ## types of DEGs
  DE_up <- round(DE * up)
  DE_down <- DE - DE_up 
  
  DE_idx <- c(rep(0, nDE), rep(1, DE_up), rep(-1, DE_down))
  DE_idx <- DE_idx[sample.int(length(DE_idx))] ## resample
  
  # h = vector indicating which pseudo-genes to re-simulate
  h <- rep(TRUE, nGenes)
  counts <- matrix(0, nrow = nGenes, ncol = 2 * m)
  
  ## log fold change, approximately half positive, half negative
  delta <- rep(0, nGenes)
  if (is.function(fc)) {
    lfc <- fc(DE)
  } else {
    lfc <- log2(fc)
  }
  delta[DE_idx != 0] <- lfc * DE_idx[DE_idx != 0]
  
  selected_genes <- true_means <- true_disps <- rep(0, nGenes)
  left_genes <- 1:length(mu)
  lambda <- phi <- matrix(0, nrow = nGenes, ncol = 2 * m)
  
  while(any(h)){
    temp <- sample.int(length(mu), sum(h), replace)
    temp <- temp[order(temp)]
    selected_genes[h] <- left_genes[temp]
    if (replace == FALSE){
      left_genes <- left_genes[-temp]
    }
    
    true_means[h] <- mu[selected_genes[h]]
    true_disps[h] <- disp[selected_genes[h]]
    
    lambda[h,] <- matrix(true_means[h], ncol = 1) %*% 
      matrix(rep(1, 2 * m), nrow = 1) * 
      cbind(matrix(rep(exp(delta[h]/2), m), ncol = m), # case
            matrix(rep(exp(-delta[h]/2), m), ncol = m)) # control
    ## mean of counts
    
    phi[h,] <- matrix(rep(true_disps[h], 2 * m), ncol = 2 * m)
    ## dispersion of counts
    
    counts[h,] <- rnegbin(sum(h) * 2 * m, lambda[h,], 1 / phi[h,])
    h <- (rowSums(cpm(counts) > 2) < 3)
    # print(sum(h))
  }
  
  if(any(rowSums(cpm(counts) > 2) < 3 ))
    print("Error: Failed to simulate data: some genes are not expressed.")
  
  list(counts = counts, 
       group = arg$group, 
       lambda0 = lambda[, 1],   # mean counts in control group
       phi0 = phi[ ,1],   # dispersion
       de = de,   # DE indicator
       delta = delta  # log2 fold change
  )
}
