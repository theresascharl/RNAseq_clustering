This repository contains R code to reproduce the analysis of the paper:
A Clustering Procedure for Three-Way RNA Sequencing Data 
Using Data Transformations and Matrix-Variate Gaussian Mixture Models

1 Figure1.R

Contains code to reproduce Figure 1:
- generation of 2 artificial data sets
- visualization in the simplex and R2

2 fission_data_preprocessing.R

Contains code for Pre-processing of RNA-seq data:

- normalised expression profiles of the genes across time points for 
  a biological unit and experiment are obtained, 
- averages are taken across biological replicates and finally 
- differentially expressed genes are identified to 
  reduce the number of observations.

Output:
- ALR profiles of differentially expressed genes
  fission_alr_de_flat.RData
- mean profiles of differentially expressed genes
  fission_mean_profiles_de_flat.RData
- ALR profiles of differentially expressed genes array
  fission_ALR_de_array.RData
  
3 Figure2.R

Contains code to reproduce Figure 2 of the preprocessing steps of 
threeway RNAseq data

4 fission_threeway_clustering.R

- select the model based on expert knowledge, for fission data it is G-VVI-VV
- select the number of clusters based on ICL

Output:
- fission_ALR_G-VVI-VV_11.RData

5 fission_twoway_clustering.R

- select the model based on expert knowledge, for fission data it is VVV
- select the number of clusters based on ICL

Output:
- fission_Malr_VVV_4.RData
- visualization of cluster solution: Figure7.pdf

6 fission_cluster_postprocessing.R

post processing of threeway clustering:
- density-based silhouette plot: Figure3.pdf
- cluster map: Figure4.pdf
- visualization of cluster solution: Figure5.pdf

post processing of twoway clustering:
- density-based silhouette plot: Figure6a.pdf
- cluster map: Figure6b.pdf
