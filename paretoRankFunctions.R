#===============================================================================
# File:    paretoRankFunctions.R
# Author:  Jeff A.Tracey
# Created: Mar 14, 2007
# Objective:
#     Pareto Ranking to Support MCDA.
#===============================================================================

library(compiler)
library(Rcpp)

dyn.load("fast_pareto_rank.so")

rowHasNoNA <- function(x) {
  res <- !any(is.na(x))
  return(res)
}

rescaleMaxToMinByCol <- function(x, maxToMin=rep(TRUE, ncol(x))) {
  x <- as.matrix(x)
  m <- ncol(x)
  res <- matrix(NA, nrow=nrow(x), ncol=m)
  for (i in 1:m) {
    if (maxToMin[i]) {
      res[,i] <- (max(x[,i], na.rm=TRUE) - x[,i])/
        (max(x[,i], na.rm=TRUE) - min(x[,i], na.rm=TRUE))
    } else {
      res[,i] <- (x[,i] - min(x[,i], na.rm=TRUE))/
        (max(x[,i], na.rm=TRUE) - min(x[,i], na.rm=TRUE))
    }
  }
  return(res)
}

cppParetoRank <- function(obj) {
  i <- which(apply(obj, 1, rowHasNoNA))
  out <- rep(NA, nrow(obj))
  res <- .Call('paretoRank', as.matrix(obj[i,]))
  out[i] <- as.vector(res)
  return(out)
}

# not really faster
cppFastParetoRank <- function(obj) {
  i <- which(apply(obj, 1, rowHasNoNA))
  out <- rep(NA, nrow(obj))
  res <- .Call('fastParetoRank', as.matrix(obj[i,]))
  out[i] <- as.vector(res)
  return(out)
}

# is really faster
cppFasterParetoRank <- function(obj) {
  i <- which(apply(obj, 1, rowHasNoNA))
  out <- rep(NA, nrow(obj))
  res <- .Call('fasterParetoRank', as.matrix(obj[i,]))
  out[i] <- as.vector(res)
  return(out)
}

# true if i dominates j
cppDominates <- function(obj, i, j, print.obj=FALSE) {
  res <- .Call('dominates', as.matrix(obj), i, j)
  if (print.obj) {
    print(obj[i,])
    print(obj[j,])
  }
  return(as.logical(res))
}

# criteria.list is a list of vectors of field names in
# order of priority
hierarchicalParetoRank <- function(obj, criteria.list) {
    res <- matrix(NA, nrow=nrow(obj), ncol=length(criteria.list))
    for (i in 1:length()) {
        res[,i] <- cppParetoRank(obj[criteria.list[[i]]])
    }
    return(res)
}

