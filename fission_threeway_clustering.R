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
for(i in 1:20) {
  print(i)
  set.seed(123)
  init <- MatTrans.init(X1, K = i, n.start = 10)
  M[[i]] <- MatTrans.EM(X1,  initial = init, model = "G-VVI-VV", 
                        row.skew = TRUE, col.skew = TRUE, 
                        trans = "None", silent = TRUE, size.control = 10)
}
save(M, file = "fission_ALR_G-VVI-VV_1to20.RData")

str(M)

# select k based on BIC ----
Mbic <- rep(0,20)
for(i in 1:20) Mbic[i] <- M[[i]]$best.bic
which.min(Mbic)
## 12 cluster according to BIC

# select k based on ILC ----
icl_M <- function(object, ...)
{
  z <- object$best.result[[1]]$gamma
  n <- nrow(z)
  if(is.null(z)) z <- matrix(1, nrow = n, ncol = 1)
  C <- matrix(0, n, ncol(z))
  for(i in 1:n) 
    C[i, which.max(z[i,])] <- 1
  object$best.bic + 2*sum(C * ifelse(z > 0, log(z), 0))
}

my_icl <- rep(NA,20)
for(i in 1:20) my_icl[i] <- try(icl_M(M[[i]])) 
which.min(as.numeric(my_icl))
## 12 cluster according to ICL

##########################################################
# three-way clustering for G-VVI-VV with 12 components ----

# errors occur in single starts and iterations but it still works
set.seed(123)
init <- MatTrans.init(X1, K = 12, n.start = 10)

M <- MatTrans.EM(X1, initial = init, model = "G-VVI-VV", 
                 row.skew = TRUE, col.skew = TRUE, 
                 trans = "None", silent = TRUE, size.control = 10)

save(M, file = "fission_ALR_G-VVI-VV_12.RData")

str(M)

M$best.bic

M$best.model

M$result[[1]]

table(M$best.result[[1]]$id)
genes[M$best.result[[1]]$id == 6]

# cluster of SPNCRNA.1164 is 6 which is cluster 10 after ordering by dbsi
M$best.result[[1]]$id[grep("SPNCRNA.1164", genes)]
cl6 <- genes[M$best.result[[1]]$id == 6]

###############################################
## visualize a subset of clusters flattened out
load("fission_ALR_G-VVI-VV_12.RData")
load("fission_alr_de_flat.RData")
load("fission_mean_profiles_de_flat.RData")
load("de_genes_fission_complete.RData")
load("mprofiles_fission.RData")

source("functions_cluster_postprocessing.R")

## get cluster info to reorder the clusters ----
cinfo <- cluster_map(M = M, genes = genes, type = "3way", 
            mprofiles = mprofiles, npcol = 4)$info 
data1 <- fission_alr_de_flat[rownames(cinfo),]
data2 <- fission_mean_profiles_de_flat[rownames(cinfo),]

col <- hcl.colors(12)
palette(col)

pdf("Figure5.pdf", width = 10, height=6)
par(mfrow=c(2,4), mar=c(1,4.1,2.1,0.1))
matplot(1:10,t(data1[cinfo$cluster=="cluster 1",]),main="cluster 1", ty="l",
        ylab="alr", col=1, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
matplot(1:10,t(data1[cinfo$cluster=="cluster 2",]),main="cluster 2", ty="l",
        ylab="alr", col=2, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
matplot(1:10,t(data1[cinfo$cluster=="cluster 7",]),main="cluster 7", ty="l",
        ylab="alr", col=6, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
matplot(1:10,t(data1[cinfo$cluster=="cluster 10",]),main="cluster 10", ty="l",
        ylab="alr", col=11, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
par(mar=c(1,4.1,0.1,0.1))

matplot(1:12,t(data2[cinfo$cluster=="cluster 1",]),main="", ty="l",
        ylab="mprofiles", col=1, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
matplot(1:12,t(data2[cinfo$cluster=="cluster 2",]),main="", ty="l",
        ylab="mprofiles", col=2, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
matplot(1:12,t(data2[cinfo$cluster=="cluster 7",]),main="", ty="l",
        ylab="mprofiles", col=6, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
matplot(1:12,t(data2[cinfo$cluster=="cluster 10",]),main="", ty="l",
        ylab="mprofiles", col=11, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
dev.off()




