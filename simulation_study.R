## load packages ----
library("matrixNormal")
library("MatTransMix") # model-based three-way clustering
library("mclust") # model-based two-way clustering 
library("cluster")    # clustering algorithms
library("factoextra") # clustering algorithms & visualization
library("e1071") # adjusted Rand index
library("abind")

source("functions_cluster_postprocessing.R")

## Step 1 - generate data ----
## Generate 100 datasets using cluster mean, covarince and size
## from the fission three-way cluster solution

## load parameters for simulated data ---
load("fission_ALR_G-VVI-VV_10.RData")

M1 <- M$best.result[[1]]
mu <- M1$Mu # means
sigma <- M1$Sigma # covariance between experiments U
psi <- M1$Psi # covariance between time points V
size <- table(M1$id) # cluster size

## generate 100 artificial data sets
set.seed(123)
for (d in 1:100) {
  X1 <- lapply(seq_along(size), function(i) {
      replicate(size[i], matrixNormal::rmatnorm(M = mu[,,i], 
                                                U = sigma[,,i], 
                                                V = psi[,,i]))
  })
  X1 <- do.call("abind", X1)
  save(X1, file = paste0("data", d, ".RData"))
}

## Step 2 - three-way clustering ----
## by model-based clustering using MatTransMix

## initialize result ---

BIC_ID <- matrix(ncol = 100, nrow = sum(size))
ICL_ID <- matrix(ncol = 100, nrow = sum(size))

for (d in 1:100) {
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
  
  set.seed(123)
  init <- MatTrans.init(X1, K = mbic, n.start = 10)

  M <- MatTrans.EM(X1, initial = init, model = "G-VVI-VV", 
                   row.skew = TRUE, col.skew = TRUE, 
                   trans = "None", silent = TRUE, size.control = 10)
  BIC_ID[, d] <- M$best.result[[1]]$id
  save(BIC_ID, file = "BIC_ID_threeway.RData")

  set.seed(123)
  init <- MatTrans.init(X1, K = micl, n.start = 10)
  
  M <- MatTrans.EM(X1, initial = init, model = "G-VVI-VV", 
                   row.skew = TRUE, col.skew = TRUE, 
                   trans = "None", silent = TRUE, size.control = 10)
  ICL_ID[, d] <- M$best.result[[1]]$id
  save(ICL_ID, file = "ICL_ID_threeway.RData")  
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
result <- matrix(NA, nrow = 600, ncol = 5)
result <- as.data.frame(result)
colnames(result) <- c("dataset", "cluster.method", "nclus.method", "crand", "K")
result$dataset <- rep(1:100, 6)
result$cluster.method <- rep(c("3-way", "2-way", "kmeans"), each = 200)
result$cluster.method <- ordered(result$cluster.method, 
                                 levels = c("3-way", "2-way", "kmeans"))
result$nclus.method <- rep(c("BIC", "ICL", "BIC", "ICL", "Silhouette",
                             "Fixed.10"), each = 100)
result$nclus.method <- ordered(result$nclus.method, 
                               levels = c("BIC", "ICL", "Silhouette", "Fixed.10"))

## load three-way BIC results
load("BIC_ID_threeway.RData")
res <- vapply(seq_len(ncol(BIC_ID)), function(i)
    compareMatchedClasses(true, BIC_ID[,i])$crand,
    numeric(1))

result$crand[1:100] <- res
result$K[1:100] <- apply(BIC_ID, 2, max)

## load three-way ICL results
load("ICL_ID_threeway.RData")
res <- vapply(seq_len(ncol(ICL_ID)), function(i)
    compareMatchedClasses(true, ICL_ID[,i])$crand,
    numeric(1))

result$crand[101:200] <- res
result$K[101:200] <- apply(ICL_ID, 2, max)

## load two-way BIC results
load("BIC_ID_twoway.RData")
res <- vapply(seq_len(ncol(BIC_ID)), function(i)
    compareMatchedClasses(true, BIC_ID[,i])$crand,
    numeric(1))

result$crand[201:300] <- res
result$K[201:300] <- apply(BIC_ID, 2, max)

## load two-way ICL results
load("ICL_ID_twoway.RData")
res <- vapply(seq_len(ncol(ICL_ID)), function(i)
    compareMatchedClasses(true, ICL_ID[,i])$crand,
    numeric(1))

result$crand[301:400] <- res
result$K[301:400] <- apply(ICL_ID, 2, max)

## load kmeans Silhouette results
load("kmeans_silhouette.RData")
res <- vapply(seq_len(ncol(SIL_ID)), function(i)
    compareMatchedClasses(true, SIL_ID[,i])$crand,
    numeric(1))

result$crand[401:500] <- res
result$K[401:500] <- apply(SIL_ID, 2, max)

## load kmeans with 10 clusters
load("kmeans_ten.RData")
res <- vapply(seq_len(ncol(TEN_ID)), function(i)
    compareMatchedClasses(true, TEN_ID[,i])$crand,
    numeric(1))

result$crand[501:600] <- res
result$K[501:600] <- apply(TEN_ID, 2, max)

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
mcols <- mycols[c(1, 1, 4, 2, 2, 3)]

clusMeth <- data.frame(cluster.method = c("3-way", "2-way", "kmeans",
                                          "3-way", "2-way", "kmeans"),
                       clmethod = c("3-way", "2-way", "kmeans",
                                    "3-way", "2-way", "kmeans"),
                       nclus.method = c("BIC","ICL", "BIC", "ICL", 
                                     "Silhouette", "Fixed.10"),
                       mcols = mcols)
clusMeth$cluster.method <- ordered(clusMeth$cluster.method, 
                                   levels = c("3-way","2-way","kmeans"))
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
  facet_wrap(method ~ cluster.method) +
  theme( strip.text.x = element_blank() ) +
  theme(legend.position = "none") +
  geom_text(data = clusMeth, aes(x = 6, y = 92, label = clmethod)) 
dev.off()

  

