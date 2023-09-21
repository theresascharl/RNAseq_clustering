## access to Bioconductor workflow
## http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

# load packages ----
library("fission")
library("coseq")
library("robCompositions")
library("DESeq2")

# load the data and check the sample info ----
data("fission", package = "fission")
tab <- fission@colData
as.data.frame(tab)

# extract the raw counts ----
fd <- fission@assays$data@listData$counts
bt <- fission@assays$data@listData$biotype

## fission data with 7039 genes and 
## 36 samples in wt/treatment, 
## time 0, 12, 30, 60, 120, 180 and 
## 3 biol replicates

exper <- c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3")
time <- c("T00", "T15", "T30", "T60", "T120", "T180")

# extract single experiment count data ----
## data needs to be rearraged for data transformations
col <- c(1, 4, 7, 10, 13, 16)
single_exp <- list()
single_exp[[1]] <- fd[, col]
single_exp[[2]] <- fd[, col+1]
single_exp[[3]] <- fd[, col+2]
single_exp[[4]] <- fd[, col+18]
single_exp[[5]] <- fd[, col+19]
single_exp[[6]] <- fd[, col+20]

save(single_exp, file = "single_exp_fission.RData")

# calculate normalized expression profiles ----
profiles <- list()
for (i in 1:length(exper)) {
  sep <- single_exp[[i]]
  counts <- sep[complete.cases(sep), ]

  set.seed(123)
  run <- coseq(counts, K = 2:20, transformation = "profile", norm = "DESeq", 
               meanFilterCutoff = 50, model = "kmeans", 
               nstart = 10, iter.max = 10)
  STprofile <- tcounts(run)
  colnames(STprofile) <- time
  
  profiles[[i]] <- na.omit(STprofile)
}

save(profiles, file = "profiles_fission.RData")

# calculate mean profiles -- sorry, very inefficient and slow
mprofiles <- list()

for (i in 1:2) {
  # only use genes that are present in all 3 biol. replicates
  un <- intersect(rownames(profiles[[1+(i-1)*3]]), rownames(profiles[[2+(i-1)*3]]))
  un <- intersect(un, rownames(profiles[[3+(i-1)*3]]))
  un <- un[!is.na(un)]
  res <- matrix(nrow = length(un), ncol = 6)
  # for each gene calculate mean profile
  for (g in 1:length(un)) {
    mat <- rbind(rbind(profiles[[1+(i-1)*3]][un[g], ], 
                       profiles[[2+(i-1)*3]][un[g], ]), 
                       profiles[[3+(i-1)*3]][un[g], ])
    res[g, ] <- as.numeric(colMeans(as.data.frame(mat), na.rm = TRUE))
  }
  res <- as.data.frame(res)
  rownames(res) <- un
  colnames(res) <- time
  mprofiles[[i]] <- res
}

save(mprofiles, file = "mprofiles_fission.RData")

#########################################################
## DESeq2 differential expression analysis --------

library("DESeq2")
library("fission")
data("fission", package = "fission")

# differentially expressed between time points
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
dds1 <- DESeq(ddsTC)
resultsNames(dds1)

# for all coefficients -- also slow
de <- NULL
for (i in 2:12) {
  resLFC <- lfcShrink(dds1, coef = i, type = "apeglm")
  de1 <- rownames(resLFC[!is.na(resLFC$padj) & resLFC$padj < 0.01 & abs(resLFC$log2FoldChange) > 1, ])
  de <- c(de, de1)
}

de <- unique(de)
save(de, file = "de_fission_all_coeff.RData")

## plot single genes with DESeq2 ----

fiss <- plotCounts(ddsTC, "SPBC2F12.09c", 
                   intgroup = c("minute", "strain"), returnData = TRUE)
fiss$minute <- as.numeric(as.character(fiss$minute))
ggplot(fiss, 
       aes(x = minute, y = count, color = strain, group = strain)) + 
  geom_point() + stat_summary(fun.y = mean, geom = "line") +
  scale_y_log10()


########################################################
## compositional data analysis ----

library("robCompositions")
load("de_fission_all_coeff.RData")

# ALR transformation of mean profiles for the differentially expressed genes
alr_object <- list()
for (i in 1:2) {
  r1 <- addLR(mprofiles[[i]], ivar = 1)
  ok <- rownames(mprofiles[[i]]) %in% de
  alr_object[[i]] <- r1$x.alr[ok, ]
}

save(alr_object, file = "trafo_alr_fission.RData")

##########################################################
## set of genes differentially expressed and present in 
## all experiments

ok <- intersect(rownames(alr_object[[1]]), rownames(alr_object[[2]]))
alr_object[[1]] <- alr_object[[1]][ok, ]
alr_object[[2]] <- alr_object[[2]][ok, ]

test <- cbind(alr_object[[1]], alr_object[[2]])
ok <- complete.cases(test) 

genes <- rownames(alr_object[[1]][ok, , drop = FALSE])
save(genes, file = "de_genes_fission_complete.RData")

#########################################################
## ALR profiles of differentially expressed genes flat

fission_alr_de_flat <- test[ok, ]
save(fission_alr_de_flat, file = "fission_alr_de_flat.RData")

#########################################################
## ALR profiles of differentially expressed genes array

n <- length(genes)
p <- 2
T <- 5
X1 <- array(NA_real_, dim = c(p, T, n))

X1[1, , ] <- t(alr_object[[1]][genes, , drop = FALSE])
X1[2, , ] <- t(alr_object[[2]][genes, , drop = FALSE])
save(X1, file = "fission_ALR_de_array.RData")


#########################################################
## mean profiles of differentially expressed genes

fission_mean_profiles_de_flat <- cbind(mprofiles[[1]][genes, ],
                                       mprofiles[[2]][genes, ])
save(fission_mean_profiles_de_flat, file = "fission_mean_profiles_de_flat.RData")

