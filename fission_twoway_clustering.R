#########################################################
## 2-way clustering of transformed data

# load packages and data ----
library("mclust")
library("RColorBrewer")
load("fission_alr_de_flat.RData")
load("fission_mean_profiles_de_flat.RData") # for visualization

## experiments flattened out to 769 genes and 10 coordinates
data1 <- fission_alr_de_flat

## only allow VVV ----
set.seed(123)
malr_VVV <- Mclust(data = data1, G = 1:20, modelNames = "VVV")
plot(malr_VVV, what = "BIC")
summary(malr_VVV)
# BIC selects 5 clusters

## select the number of clusters by ICL
Mmalr <- vector("list", 20)
set.seed(123)
for (i in 1:20) {
    Mmalr[[i]] <- Mclust(data = data1, G = i, modelNames = "VVV")
}
# no results for more than 14 clusters
Micl <- vapply(Mmalr, function(x) tryCatch(icl(x), error = function(e) NA_real_), 
               numeric(1))
which.max(Micl)
# ICL would choose 4 clusters for ALR as we are searching for the max here due
# to the mclust implementation

set.seed(123)
malr_VVV <- Mclust(data = data1, G = 4, modelNames = "VVV")
save(malr_VVV, file = "fission_Malr_VVV_4.RData")

