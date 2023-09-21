## load packages ----
library("tidyverse")
library("robCompositions")
library("flexclust")
library("RColorBrewer")
library("ggpubr")

## source helper functions ----
source("functions_cluster_postprocessing.R")

## select as many colors as clusters ------------------
col <- hcl.colors(12)
palette(col)

## load threeway cluster result ----
load("fission_ALR_G-VVI-VV_12.RData")
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


###########################################################
## compare to two-way clustering ----

load("fission_Malr_VVV_4.RData")

## calculate the adjusted rand index ----
library("e1071")
compareMatchedClasses(malr_VVV$classification, M$best.result[[1]]$id)$crand

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

