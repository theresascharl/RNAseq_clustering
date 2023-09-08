#########################################################
## 2-way clustering of transformed data

# load packages and data ----
library("mclust")
library("RColorBrewer")
load("fission_alr_de_flat.RData")
load("fission_mean_profiles_de_flat.RData") # for visualization

## experiments flatten out to 772 genes and 10 coordinates
data1 <- fission_alr_de_flat

## only allow VVV ----
set.seed(123)
malr_VVV <- Mclust(data = data1, G = 1:20, modelNames = "VVV")
mclust::plot.Mclust(malr_VVV, what="BIC")
summary(malr_VVV)
# BIC selects 5 clusters

Mmalr <- list()
set.seed(123)
for(i in 1:20) Mmalr[[i]] <- Mclust(data = data1, G = i, modelNames = "VVV")

Micl <- rep(NA,20)
for(i in 1:20) Micl[i] <- icl(Mmalr[[i]])
which.mac(Micl)
# ICL would choose 4 clusters for ALR

set.seed(123)
malr_VVV <- Mclust(data = data1, G = 4, modelNames = "VVV")
save(malr_VVV, file = "fission_Malr_VVV_4.RData")

###########################################
## visualization of 2D cluster solution----

# colors
col <- hcl.colors(4)
palette(col)

# original data flattened out
data2 <- fission_mean_profiles_de_flat

pdf("Figure7.pdf", width = 10, height=6)
par(mfrow=c(2,4), mar=c(1,4.1,0.1,0.1))
matplot(1:10,t(data1[malr_VVV$classification==3,]),main="", ty="l",
        ylab="alr", col=1, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
matplot(1:10,t(data1[malr_VVV$classification==4,]),main="", ty="l",
        ylab="alr", col=2, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
matplot(1:10,t(data1[malr_VVV$classification==2,]),main="", ty="l",
        ylab="alr", col=3, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))
matplot(1:10,t(data1[malr_VVV$classification==1,]),main="", ty="l",
        ylab="alr", col=4, lty=1, xlab="",ylim=c(-3,5.5), xaxt="n")
abline(v=c(5.5))

matplot(1:12,t(data2[malr_VVV$classification==3,]),main="", ty="l",
        ylab="mprofiles", col=1, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
matplot(1:12,t(data2[malr_VVV$classification==4,]),main="", ty="l",
        ylab="mprofiles", col=2, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
matplot(1:12,t(data2[malr_VVV$classification==2,]),main="", ty="l",
        ylab="mprofiles", col=3, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
matplot(1:12,t(data2[malr_VVV$classification==1,]),main="", ty="l",
        ylab="mprofiles", col=4, lty=1, xlab="",ylim=c(0,1), xaxt="n")
abline(v=c(6.5))
dev.off()
