This repository contains R code to reproduce the analysis of the
paper:

A Clustering Procedure for Three-Way RNA Sequencing Data Using Data
Transformations and Matrix-Variate Gaussian Mixture Models

1 Figure1.R

Contains code to reproduce Figure 1:
- Generate the 2 artificial data sets,
- Visualise the data in the simplex and R2.

2 fission_data_preprocessing.R

Contains code for Step 1: Pre-processing of RNA-seq data:

- Obtain normalised expression profiles of the genes across time
  points for a biological unit and experiment.
- Take averages across biological replicates.
- Identify differentially expressed genes to reduce the number of
  observations.

Output:
- ALR profiles of differentially expressed genes
  fission_alr_de_flat.RData
- mean profiles of differentially expressed genes
  fission_mean_profiles_de_flat.RData
- ALR profiles of differentially expressed genes array
  fission_ALR_de_array.RData
  
3 Figure2.R

Contains code to reproduce Figure 2 of the pre-processing steps of
three-way RNAseq data.

4 fission_threeway_clustering.R

- Fit finite mixture models with matrix-normal components where the
  parameters are specified as G-VVI-VV and different number of
  components.
- Select the number of clusters based on ICL.

Output:
- fission_ALR_G-VVI-VV_11.RData

5 fission_twoway_clustering.R

- Fit finite mixture models with multivariate normal components where
  the variance-covariance matrices are specified as VVV and different
  number of components.
- Select the number of clusters based on ICL.

Output:
- fission_Malr_VVV_4.RData
- visualization of cluster solution: Figure7.pdf

6 fission_cluster_postprocessing.R

Contains code for post-processing of the three-way clustering:
- dbsi information plot: Figure3.pdf
- cluster map: Figure4.pdf
- visualisation of cluster solution: Figure5.pdf

Contains code for post-processing of the two-way clustering:
- dbsi information plot: Figure6a.pdf
- cluster map: Figure6b.pdf
