## load packages ----
require("tidyverse")
require("robCompositions")
require("flexclust")
require("RColorBrewer")
require("ggpubr")

## source helper functions ----
source("functions_cluster_postprocessing.R")

## select as many colors as clusters ------------------
col <- hcl.colors(11)
col1 <- alpha(col, alpha = 0.3)
palette(col)

## load cluster result ----
load("fission_ALR_G-VVI-VV_11.RData")
load("de_genes_fission_complete.RData")
load("mprofiles_fission.RData")

## calculate dbsi ----
res <- get.si(M, genes = genes)

## silhouette plot ----
pdf("Figure3.pdf", width = 6, height = 6)
silplot_pos(res$result)
dev.off()

## cluster map ----
pdf("Figure4.pdf", width = 6, height = 6)
cluster_map(M = M, genes = genes, type = "3way", 
            mprofiles = mprofiles, npcol = 4) 
dev.off()


###########################################################
## compare to two-way clustering ----

load("fission_Malr_VVV_4.RData")

## calculate the adjusted rand index ----
library("e1071")
compareMatchedClasses(malr_VVV$classification, M$best.result[[1]]$id)$crand

## get dbsi for mclust ----
res <- get.si(malr_VVV, genes = genes, type = "2way")

pdf("Figure6a.pdf", width = 6, height = 6)
silplot_pos(res$result)
dev.off()

## cluster map ----
pdf("Figure6b.pdf", width = 6, height = 6)
cluster_map(M = malr_VVV, genes = genes, type = "2way", 
            mprofiles = mprofiles, npcol = 2) 
dev.off()

