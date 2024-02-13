## load packages ----
library("matrixNormal")
library("MatTransMix") # model-based three-way clustering
library("mclust") # model-based two-way clustering 
library("cluster")    # clustering algorithms
library("factoextra") # clustering algorithms & visualization
library("e1071") # adjusted Rand index
library("abind")

source("functions_cluster_postprocessing.R")
source("Mat.EM.R")

## Step 1 - generate data ----
## Generate 100 datasets using cluster mean, covariance and size
## from the fission three-way cluster solution with AR1 column-wise
## covariance matrix

## load parameters for simulated data ---
load("fission_ALR_G-VVI-VAR1_1to20.RData")

M1 <- fit_ar1[[10]]
mu <- M1$Mu # means
sigma <- M1$Sigma # covariance between experiments U
psi <- M1$Psi # covariance between time points V
size <- table(M1$id) # cluster size

## generate 100 artificial data sets
set.seed(123)
for (d in 1:100) {
  FILE <- paste0("data", d, ".RData")
  if (!file.exists(FILE)) {
    X1 <- lapply(seq_along(size), function(i) {
      replicate(size[i], matrixNormal::rmatnorm(M = mu[,,i], 
                                                U = sigma[,,i], 
                                                V = psi[,,i]))
    })
    X1 <- do.call("abind", X1)
    save(X1, file = FILE)
  }
}

## Step 2 - three-way clustering ----
## by model-based clustering using MatTransMix

## initialize result ---

BIC_ID <- matrix(ncol = 100, nrow = sum(size))
ICL_ID <- matrix(ncol = 100, nrow = sum(size))
BIC_AR1_ID <- matrix(ncol = 100, nrow = sum(size))
ICL_AR1_ID <- matrix(ncol = 100, nrow = sum(size))

FILE_BIC <- "BIC_ID_threeway.RData"
FILE_ICL <- "ICL_ID_threeway.RData"
FILE_AR1_BIC <- "BIC_AR1_ID_threeway.RData"
FILE_AR1_ICL <- "ICL_AR1_ID_threeway.RData"

runs <- 1:100
if (file.exists(FILE_BIC) & file.exists(FILE_ICL) &
    file.exists(FILE_AR1_BIC) & file.exists(FILE_AR1_ICL) ) {
  load(FILE_BIC)
  load(FILE_ICL)
  load(FILE_AR1_BIC)
  load(FILE_AR1_ICL)
  already_done_bic <- which(colSums(is.na(BIC_ID)) == 0)
  already_done_icl <- which(colSums(is.na(ICL_ID)) == 0)
  already_done <- intersect(already_done_bic, already_done_icl)
  runs <- setdiff(runs, already_done)
}

for (d in runs) {
  cat(paste("Dataset", d, "\n"))
  load(paste0("data", d, ".RData"))
  
  M <- vector("list", 20)
  for (i in 1:20) {
    cat("i =", i, "\n")
    set.seed(123)
    init <- MatTrans.init(X1, K = i, n.start = 10)
    M[[i]] <- MatTrans.EM(X1, initial = init, model = "G-VVI-VV", 
                          row.skew = TRUE, col.skew = TRUE, 
                          trans = "None", silent = TRUE, size.control = 10)
  }
  
  # select k based on BIC ---
  Mbic <- vapply(M, function(x) x$best.bic, numeric(1))
  mbic <- which.min(Mbic)
  cat("For dataset", d, "the minimum BIC is", mbic, ".")
  # select k based on ILC ---
  my_icl <- vapply(M, function(x) tryCatch(icl_M(x), error = function(e) NA_real_),
                   numeric(1))
  micl <- which.min(as.numeric(my_icl))
  cat("For dataset", d, "the minimum ICL is", micl, ".")
  
  BIC_ID[, d] <- M[[mbic]]$best.result[[1]]$id
  save(BIC_ID, file = FILE_BIC)
  
  ICL_ID[, d] <- M[[micl]]$best.result[[1]]$id
  save(ICL_ID, file = FILE_ICL)  
  
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
  
  BIC <- data.frame(model = "AR1", G = 1:20, BIC = vapply(fit_ar1, function(x) 
    ifelse(is.null(x$bic), NA_real_, x$bic), numeric(1)))
  
  # no results for more than 14 components
  BIC_AR1_ID[, d] <- fit_ar1[[which.min(BIC[1:14,3])]]$id
  save(BIC_AR1_ID, file = FILE_AR1_BIC)
  
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
  
  ICL_AR1_ID[, d] <- fit_ar1[[which.min(ICL[1:14,3])]]$id
  save(ICL_AR1_ID, file = FILE_AR1_ICL)
  
}


## Step 3 - two-way clustering ----
## by model-based clustering using mclust

## initialize result ---

BIC_ID <- matrix(ncol = 100, nrow = 769)
ICL_ID <- matrix(ncol = 100, nrow = 769)

for (d in 1:100) {
  cat(paste("Dataset", d, "\n"))
  load(paste0("data", d, ".RData"))
  
  data1 <- t(rbind(X1[1,,], X1[2,,]))
  
  ## only allow VVV ---
  set.seed(123)
  m_VVV <- Mclust(data = data1, G = 1:20, modelNames = "VVV")
  kbic <- m_VVV$G
  cat("For dataset", d, "the minimum BIC is", kbic, ".")

  ## select the number of clusters by ICL
  Mm <- vector("list", 20)
  set.seed(123)
  for (i in 1:20) {
    Mm[[i]] <- Mclust(data = data1, G = i, modelNames = "VVV")
  }
  # no results for more than 14 clusters
  Micl <- vapply(Mm, function(x) tryCatch(icl(x), error = function(e) NA_real_),
                 numeric(1))
  kicl <- which.max(Micl)
  cat("For dataset", d, "the minimum ICL is", kicl, ".")

  set.seed(123)
  m_VVV <- Mclust(data = data1, G = kbic, modelNames = "VVV")

  BIC_ID[, d] <- m_VVV$classification
  save(BIC_ID, file = "BIC_ID_twoway.RData")

  set.seed(123)
  m_VVV <- Mclust(data = data1, G = kicl, modelNames = "VVV")

  ICL_ID[, d] <- m_VVV$classification
  save(ICL_ID, file = "ICL_ID_twoway.RData")
}

## Step 4 - kmeans ----
## select k using average silhouette method 
## and the true number of clusters 10

## initialize result ---

SIL_ID <- matrix(ncol = 100, nrow = 769)
TEN_ID <- matrix(ncol = 100, nrow = 769)

for (z in 1:100) {
  cat(paste("Dataset", z, "\n"))
  load(paste0("data", z, ".RData"))
  
  data1 <- t(rbind(X1[1,,], X1[2,,]))
  
  ## find k using silhouette method
  
  set.seed(123)
  res <- fviz_nbclust(data1, kmeans, k.max = 20, 
                      method = "silhouette")
  silk <- which.max(res$data$y)
  
  set.seed(123)
  km_sil <- kmeans(data1, centers=silk, iter.max = 100, nstart = 10)
  print(table(km_sil$cluster))  
  
  ## use 10 clusters 
  set.seed(123)
  km_10 <- kmeans(data1, centers = 10, iter.max = 100, nstart = 10)
  
  SIL_ID[, z] <- km_sil$cluster
  save(SIL_ID, file = "kmeans_silhouette.RData")
  
  TEN_ID[, z] <- km_10$cluster
  save(TEN_ID, file = "kmeans_ten.RData")
  
}

## Step 5 - Summary ----

## true cluster membership by design
true <- rep(1:length(size), size)

## initialize
result <- matrix(NA_real_, nrow = 800, ncol = 5)
result <- as.data.frame(result)
colnames(result) <- c("dataset", "cluster.method", "nclus.method", "crand", "K")
result$dataset <- rep(1:100, 8)
result$cluster.method <- rep(c("3-way", "3-way AR1", "2-way", "kmeans"), 
                             each = 200)
result$cluster.method <- ordered(result$cluster.method, 
                                 levels = c("3-way", "3-way AR1", "2-way",
                                            "kmeans"))
result$nclus.method <- rep(c("BIC", "ICL", "BIC", "ICL", "BIC", "ICL", 
                             "Silhouette", "Fixed.10"), each = 100)
result$nclus.method <- ordered(result$nclus.method, 
                               levels = c("BIC", "ICL", "Silhouette", "Fixed.10"))

## load three-way BIC results MatTransMix
load("BIC_ID_threeway.RData")
res <- vapply(seq_len(ncol(BIC_ID)), function(i)
  compareMatchedClasses(true, BIC_ID[,i])$crand,
  numeric(1))

result$crand[1:100] <- res
result$K[1:100] <- apply(BIC_ID, 2, max)

## load three-way ICL results MatTransMix
load("ICL_ID_threeway.RData")
res <- vapply(seq_len(ncol(ICL_ID)), function(i)
  compareMatchedClasses(true, ICL_ID[,i])$crand,
  numeric(1))

result$crand[101:200] <- res
result$K[101:200] <- apply(ICL_ID, 2, max)

## load three-way BIC results Mat.EM AR1
load("BIC_AR1_ID_threeway.RData")
res <- vapply(seq_len(ncol(BIC_AR1_ID)), function(i)
  compareMatchedClasses(true, BIC_AR1_ID[,i])$crand,
  numeric(1))

result$crand[201:300] <- res
result$K[201:300] <- apply(BIC_AR1_ID, 2, max)

## load three-way ICL results Mat.EM AR1
load("ICL_AR1_ID_threeway.RData")
res <- vapply(seq_len(ncol(ICL_AR1_ID)), function(i)
  compareMatchedClasses(true, ICL_AR1_ID[,i])$crand,
  numeric(1))

result$crand[301:400] <- res
result$K[301:400] <- apply(ICL_AR1_ID, 2, max)

## load two-way BIC results
load("BIC_ID_twoway.RData")
res <- vapply(seq_len(ncol(BIC_ID)), function(i)
  compareMatchedClasses(true, BIC_ID[,i])$crand,
  numeric(1))

result$crand[401:500] <- res
result$K[401:500] <- apply(BIC_ID, 2, max)

## load two-way ICL results
load("ICL_ID_twoway.RData")
res <- vapply(seq_len(ncol(ICL_ID)), function(i)
  compareMatchedClasses(true, ICL_ID[,i])$crand,
  numeric(1))

result$crand[501:600] <- res
result$K[501:600] <- apply(ICL_ID, 2, max)

## load kmeans Silhouette results
load("kmeans_silhouette.RData")
res <- vapply(seq_len(ncol(SIL_ID)), function(i)
  compareMatchedClasses(true, SIL_ID[,i])$crand,
  numeric(1))

result$crand[601:700] <- res
result$K[601:700] <- apply(SIL_ID, 2, max)

## load kmeans with 10 clusters
load("kmeans_ten.RData")
res <- vapply(seq_len(ncol(TEN_ID)), function(i)
  compareMatchedClasses(true, TEN_ID[,i])$crand,
  numeric(1))

result$crand[701:800] <- res
result$K[701:800] <- apply(TEN_ID, 2, max)

result$method <- rep(c("A", "B"), each = 100)

save(result, file = "simulation_result.RData")

## Step 6 - visualization ----
## Visualization
library("ggplot2")
library("wesanderson")

load("simulation_result.RData")

pdf("Figure8a.pdf")
ggplot(result, aes(x = cluster.method, y = crand, fill = nclus.method)) + 
  theme_bw() + ylim(0:1) +
  scale_fill_manual(values = wes_palette(n = 4, name = "GrandBudapest1")) +
  xlab("cluster method") +
  ylab("adjusted Rand index") +
  labs(fill = "Select K by") +
  geom_boxplot()
dev.off()

mycols <- wes_palette(n = 4, name = "GrandBudapest1")
mcols <- mycols[c(1, 1, 1, 4, 2, 2, 2, 3)]

clusMeth <- data.frame(cluster.method = c("3-way", "3-way AR1", "2-way", "kmeans",
                                          "3-way", "3-way AR1", "2-way", "kmeans"),
                       clmethod = c("3-way", "3-way AR1", "2-way", "kmeans",
                                    "3-way", "3-way AR1", "2-way", "kmeans"),
                       nclus.method = c("BIC","ICL", "BIC", "ICL", "BIC", "ICL",
                                        "Silhouette", "Fixed.10"),
                       mcols = mcols)
clusMeth$cluster.method <- ordered(clusMeth$cluster.method, 
                                   levels = c("3-way", "3-way AR1", "2-way",
                                              "kmeans"))
clusMeth$nclus.method <- ordered(clusMeth$nclus.method, 
                                 levels = c("BIC","ICL", "Silhouette",
                                            "Fixed.10"))
pdf("Figure8b.pdf")
ggplot(result, aes(x = K, fill = nclus.method)) + 
  geom_bar() +
  scale_fill_manual(values = wes_palette(n = 4, name = "GrandBudapest1")) +
  theme_bw() +
  labs(fill = "Select K by") +
  xlab("Number of clusters K") +
  ylab("count") +
  facet_wrap(method ~ cluster.method, nrow = 2) +
  theme( strip.text.x = element_blank() ) +
  theme(legend.position = "none") +
  geom_text(data = clusMeth, aes(x = 6, y = 92, label = clmethod)) 
dev.off()



