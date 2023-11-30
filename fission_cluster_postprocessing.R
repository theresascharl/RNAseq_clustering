## load packages ----
library("tidyverse")
library("robCompositions")
library("flexclust")
library("RColorBrewer")
library("ggpubr")

## source helper functions ----
source("functions_cluster_postprocessing.R")

## select as many colors as clusters ------------------
col <- hcl.colors(10)
palette(col)

## load threeway cluster result ----
load("fission_ALR_G-VVI-VV_10.RData")
load("de_genes_fission_complete.RData")
load("mprofiles_fission.RData")

## calculate dbsi ----
dbsi_3way <- get_dbsi(M, genes = genes)

## dbsi plot ----
pdf("Figure3.pdf", width = 6, height = 6)
dbsi_plot(dbsi_3way$result)
dev.off()

## cluster map ----
pdf("Figure4.pdf", width = 6, height = 6)
cluster_map(M = M, genes = genes, type = "3way", 
            mprofiles = mprofiles, npcol = 4)$plot 
dev.off()

###############################################
## visualize a subset of clusters flattened out
load("fission_ALR_G-VVI-VV_10.RData")
load("fission_alr_de_flat.RData")
load("fission_mean_profiles_de_flat.RData")
load("de_genes_fission_complete.RData")
load("mprofiles_fission.RData")

source("functions_cluster_postprocessing.R")

## get cluster info to reorder the clusters ----
cinfo <- cluster_map(M = M, genes = genes, type = "3way", 
                     mprofiles = mprofiles, npcol = 4)$info 
data1 <- fission_alr_de_flat[rownames(cinfo), ]
data2 <- fission_mean_profiles_de_flat[rownames(cinfo), ]

col <- hcl.colors(10)
palette(col)

pdf("Figure5.pdf", width = 10, height = 6)
par(mfrow = c(2, 4), mar = c(3.1, 4.1, 2.1, 0.1))
matplot(1:10, t(data1[cinfo$cluster == "cluster 1", ]), main = "cluster 1", 
        ty = "n", ylab = "alr", col = 1, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[cinfo$cluster == "cluster 1", ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 1)
  lines(6:10, dat[i,6:10], col = 1)
}
abline(v = 5.5)
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)

matplot(1:10, t(data1[cinfo$cluster == "cluster 2", ]), main = "cluster 2", 
        ty = "n", ylab = "alr", col = 2, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[cinfo$cluster == "cluster 2", ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 2)
  lines(6:10, dat[i,6:10], col = 2)
}
abline(v = 5.5)
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)

matplot(1:10, t(data1[cinfo$cluster == "cluster 5", ]), main = "cluster 5", 
        ty = "n", ylab = "alr", col = 5, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[cinfo$cluster == "cluster 5", ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 5)
  lines(6:10, dat[i,6:10], col = 5)
}
abline(v = 5.5)
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)

matplot(1:10, t(data1[cinfo$cluster == "cluster 8", ]), main = "cluster 8", 
        ty = "n", ylab = "alr", col = 8, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[cinfo$cluster == "cluster 8", ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 8)
  lines(6:10, dat[i,6:10], col = 8)
}
abline(v = 5.5)
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)
par(mar = c(1, 4.1, 0.1, 0.1))

matplot(1:12, t(data2[cinfo$cluster == "cluster 1", ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 1, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[cinfo$cluster == "cluster 1", ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 1)
  lines(7:12, dat[i,7:12], col = 1)
}
abline(v = 6.5)

matplot(1:12, t(data2[cinfo$cluster == "cluster 2", ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 2, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[cinfo$cluster == "cluster 2", ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 2)
  lines(7:12, dat[i,7:12], col = 2)
}
abline(v = 6.5)

matplot(1:12, t(data2[cinfo$cluster == "cluster 5", ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 5, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[cinfo$cluster == "cluster 5", ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 5)
  lines(7:12, dat[i,7:12], col = 5)
}
abline(v = 6.5)

matplot(1:12, t(data2[cinfo$cluster == "cluster 8", ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 8, lty = 1, xlab = "", ylim = c(0, 1), xaxt = "n")
dat <- data2[cinfo$cluster == "cluster 8", ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 8)
  lines(7:12, dat[i,7:12], col = 8)
}
abline(v = 6.5)
dev.off()

###########################################################
## compare to two-way clustering ----

load("fission_Malr_VVV_4.RData")

## calculate the adjusted rand index ----
library("e1071")
compareMatchedClasses(M$best.result[[1]]$id, malr_VVV$classification)$crand
# 0.15 compare large to small
compareMatchedClasses(malr_VVV$classification, M$best.result[[1]]$id)$crand
# 0.78 compare small to large 
(tb <- table(threeway = M$best.result[[1]]$id, twoway = malr_VVV$classification))

library("xtable")
xtable(tb)

## now also compare to kmeans
load("fission_kmeans_3_silhouette.RData")

compareMatchedClasses(M$best.result[[1]]$id, kmeans3)$crand
#0.14
(tb <- table(threeway = M$best.result[[1]]$id, kmeans = kmeans3))
xtable(tb)

## get dbsi for mclust ----
dbsi_2way <- get_dbsi(malr_VVV, genes = genes, type = "2way")

pdf("Figure6a.pdf", width = 6, height = 6)
dbsi_plot(dbsi_2way$result)
dev.off()

## cluster map ----
pdf("Figure6b.pdf", width = 6, height = 6)
cluster_map(M = malr_VVV, genes = genes, type = "2way", 
            mprofiles = mprofiles, npcol = 2)$plot 
dev.off()

###########################################
## visualization of 2D cluster solution----
load("fission_alr_de_flat.RData")
load("fission_mean_profiles_de_flat.RData") # for visualization

## experiments flattened out to 769 genes and 10 coordinates
data1 <- fission_alr_de_flat

# colors
col <- hcl.colors(4)
palette(col)

# original data flattened out
data2 <- fission_mean_profiles_de_flat

# clusters need to be rearranged to match the order by average dbsi
pdf("Figure7.pdf", width = 10, height = 6)
par(mfrow = c(2, 4), mar = c(3.1, 4.1, 2.1, 0.1))
matplot(1:10, t(data1[malr_VVV$classification == 1, ]), main = "cluster 1", 
        ty = "n", ylab = "alr", col = 1, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[malr_VVV$classification == 1, ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 1)
  lines(6:10, dat[i,6:10], col = 1)
}
abline(v = c(5.5))
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)
matplot(1:10, t(data1[malr_VVV$classification == 4, ]), main = "cluster 2", 
        ty = "n", ylab = "alr", col = 2, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[malr_VVV$classification == 4, ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 2)
  lines(6:10, dat[i,6:10], col = 2)
}
abline(v = c(5.5))
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)
matplot(1:10, t(data1[malr_VVV$classification == 2, ]), main = "cluster 3", 
        ty = "n", ylab = "alr", col = 3, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[malr_VVV$classification == 2, ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 3)
  lines(6:10, dat[i,6:10], col = 3)
}
abline(v = c(5.5))
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)
matplot(1:10, t(data1[malr_VVV$classification == 3, ]), main = "cluster 4", 
        ty = "n", ylab = "alr", col = 4, lty = 1, xlab = "", ylim = c(-3, 5.5),
        xaxt = "n")
dat <- data1[malr_VVV$classification == 3, ]
for(i in 1:nrow(dat)) {
  lines(1:5, dat[i,1:5], col = 4)
  lines(6:10, dat[i,6:10], col = 4)
}
abline(v = c(5.5))
axis(1, at = c(3,8), label = c("WT", "Mut"), tick = FALSE)
par(mar = c(1, 4.1, 0.1, 0.1))

matplot(1:12, t(data2[malr_VVV$classification == 1, ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 1, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[malr_VVV$classification == 1, ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 1)
  lines(7:12, dat[i,7:12], col = 1)
}
abline(v = c(6.5))
matplot(1:12, t(data2[malr_VVV$classification == 4, ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 2, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[malr_VVV$classification == 4, ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 2)
  lines(7:12, dat[i,7:12], col = 2)
}
abline(v = c(6.5))
matplot(1:12, t(data2[malr_VVV$classification == 2, ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 3, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[malr_VVV$classification == 2, ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 3)
  lines(7:12, dat[i,7:12], col = 3)
}
abline(v = c(6.5))
matplot(1:12, t(data2[malr_VVV$classification == 3, ]), main = "", ty = "n", 
        ylab = "mprofiles", col = 4, lty = 1, xlab = "", ylim = c(0, 1), 
        xaxt = "n")
dat <- data2[malr_VVV$classification == 3, ]
for(i in 1:nrow(dat)) {
  lines(1:6, dat[i,1:6], col = 4)
  lines(7:12, dat[i,7:12], col = 4)
}
abline(v = c(6.5))
dev.off()

