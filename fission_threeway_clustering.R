#########################################################
## 3-way clustering

## ALR

library("MatTransMix")
load("fission_ALR_de_array.RData") ## load array X1
load("de_genes_fission_complete.RData") ## load gene vector

## data arrangement for 3-way clustering
# n number of genes; p number of experiments; T number of coordinates, times

n <- length(genes)
p <- 2
T <- 5
str(X1)

#########################################################
## models for different numbers of clusters

# use G-VVI-VV as it represents the structure of the experiment
# find the optimal number of clusters for G-VVI-VV ----

M <- list()
for (i in 1:20) {
  cat("i =", i, "\n")
  set.seed(123)
  init <- MatTrans.init(X1, K = i, n.start = 10)
  M[[i]] <- MatTrans.EM(X1, initial = init, model = "G-VVI-VV", 
                        row.skew = TRUE, col.skew = TRUE, 
                        trans = "None", silent = TRUE, size.control = 10)
}
save(M, file = "fission_ALR_G-VVI-VV_1to20.RData")


#########################################################
## refinement step where the structure of the column-wise covariance
## matrix is restricted to AR1

set.seed(1234)
require("MatTransMix")
source("Mat.EM.R")

load("fission_ALR_G-VVI-VV_1to20.RData")
load("fission_ALR_de_array.RData")

fit_ar1 <- vector("list", 20)
control <- list(verbose = 1000, iter.max = 1000)
  
for (g in 1:20) {
    cat("g = ", g, "\n")
    if (length(M[[g]]$result) > 0) {
        initial <- M[[g]]$result[[1]]
        fit_ar1[[g]] <- Mat.EM(X1, initial,
                               row_model = "diag",
                               col_model = "ar1",
                               control = control)
    }
}
save(fit_ar1, file = "fission_ALR_G-VVI-VAR1_1to20.RData")

## get BIC of MatTransMix object M with full column-wise covariance matrix
bic <- sapply(M, "[[", "best.bic")
bic <- ifelse(is.finite(bic), bic, NA_real_)

## get BIC of refinement step
BIC <- data.frame(model = "AR1", G = 1:20, BIC = vapply(fit_ar1, function(x) 
  ifelse(is.null(x$bic), NA_real_, x$bic), numeric(1)))

which.min(BIC[1:14, 3])
min(BIC[1:14, 3]) ## AR1 gives smaller BIC than original model

# select k based on ILC ----
icl_M <- function(object, ...)
{
  z <- object$best.result[[1]]$gamma
  n <- nrow(z)
  if (is.null(z)) {
    z <- matrix(1, nrow = n, ncol = 1)
  }
  C <- matrix(0, n, ncol(z))
  for (i in 1:n) {
    C[i, which.max(z[i,])] <- 1
  }
  object$best.bic + 2 * sum(C * ifelse(z > 0, log(z), 0))
}

icl <- vapply(M, function(x) tryCatch(icl_M(x), error = function(e) NA_real_),
              numeric(1))
which.min(as.numeric(icl))
## 10 cluster according to ICL on MatTransMix object

getICL <- function(x) {
  if (is.null(x$bic)) {
    return(NA_real_)
  }
  z <- x$gamma
  n <- nrow(z)
  if (is.null(z)) {
    z <- matrix(1, nrow = n, ncol = 1)
  }
  C <- matrix(0, n, ncol(z))
  for (i in 1:n) {
    C[i, which.max(z[i,])] <- 1
  }
  x$bic + 2 * sum(C * ifelse(z > 0, log(z), 0))
}

ICL <- data.frame(model = "AR1", G = 1:20,
                  ICL = vapply(fit_ar1, getICL, numeric(1)))
which.min(ICL[1:14, 3])
# 14 clusters according to ICL for AR1 model

