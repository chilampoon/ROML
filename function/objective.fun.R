# Optimize final probabilities using objective function for model 3
opti.RFs <- function(prs.table, grp.names) {
  suppressPackageStartupMessages(library(nloptr))
  
  ngrp <- length(grp.names)
  lb <- rep(0, ngrp) # lower bound
  ub <- rep(1, ngrp) # upper bound
  x0 <- setNames(rep(1/ngrp, ngrp), grp.names) # vector with starting values for the optimization
  
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
  
  # Rewrite the equality as 1 - (P1 + P2 + P3 + P4 + P5) = 0
  eq_c <- function(x, pw.p.mat) 1 - sum(x) 
  
  # Rewrite the inequality as P1 * P2 * P3 * P4 * P5 - 1 <= 0
  ineq_c <- function(x, pw.p.mat) prod(x) - 1
  
  # Function to transform a row in prs.table to a matrix
  rowToMat <- function(row, grp.names) {
    row.mat <- setNames(data.frame(matrix(ncol=ngrp, nrow=ngrp), stringsAsFactors=F), grp.names)
    rownames(row.mat) <- colnames(row.mat)
    for (g in 1:(ngrp-1)) {
      grp1 <- grp.names[g]
      for (r in (g+1):ngrp) {
        grp2 <- grp.names[r]
        grp <- sprintf('%s_%s', grp1, grp2)
        row.mat[g, r] <- row[, sprintf('%s.%s', grp1, grp)]
        row.mat[r, g] <- row[, sprintf('%s.%s', grp2, grp)]
      }
    }
    row.mat
  }
  
  
  # Summarize final probabilities
  final.prs <- setNames(data.frame(matrix(ncol = ngrp, nrow = 0), stringsAsFactors=F), grp.names)
  for (n in 1:nrow(prs.table)) {
    pw.p.mat <- rowToMat(prs.table[n, ], grp.names)
    opt <- nloptr(x0=x0, eval_f=obj_f, lb=lb, ub=ub, eval_g_ineq=ineq_c, eval_g_eq=eq_c, 
                  opts=list(algorithm="NLOPT_GN_ISRES", "maxeval" = 150000), pw.p.mat=pw.p.mat)
    final.prs[nrow(final.prs)+1, ] <- opt$solution
  }
  rownames(final.prs) <- rownames(prs.table)
  final.prs
}

