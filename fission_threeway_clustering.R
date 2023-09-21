#########################################################
## 3-way clustering

## ALR

library("MatTransMix")
load("fission_ALR_de_array.RData")
load("de_genes_fission_complete.RData")

## data arrangement for 3-way clustering
# n number of genes; p number of experiments; T number of coordinates, times

n <- length(genes)
p <- 2
T <- 5
str(X1)

###########################################
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

str(M)

# select k based on BIC ----
Mbic <- vapply(M, function(x) x$best.bic, numeric(1))
which.min(Mbic)
## 10 cluster according to BIC

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

my_icl <- vapply(M, function(x) tryCatch(icl_M(x), error = function(e) NA_real_),
                 numeric(1))
which.min(as.numeric(my_icl))
## 10 cluster according to ICL

##########################################################
# three-way clustering for G-VVI-VV with 10 components ----

# errors occur in single starts and iterations but it still works
set.seed(123)
init <- MatTrans.init(X1, K = 10, n.start = 10)

M <- MatTrans.EM(X1, initial = init, model = "G-VVI-VV", 
                 row.skew = TRUE, col.skew = TRUE, 
                 trans = "None", silent = TRUE, size.control = 10)

save(M, file = "fission_ALR_G-VVI-VV_10.RData")

str(M)

M$best.bic

M$best.model

M$result[[1]]

table(M$best.result[[1]]$id)

# cluster of SPNCRNA.1164 is 10 which is cluster 8 after ordering by dbsi
M$best.result[[1]]$id[grep("SPNCRNA.1164", genes)]
cl10 <- genes[M$best.result[[1]]$id == 10]
cl10

