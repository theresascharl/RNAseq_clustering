library("coseq")

###############################################
## visualize raw and transformed data

load("single_exp_fission.RData")
load("profiles_fission.RData")
load("trafo_alr_fission.RData")
load("mprofiles_fission.RData")

################################################

library("RColorBrewer")
cols <- brewer.pal(5, "Dark2")

#' Plot raw counts or profiles of single experiments
#'
#' @param data A list of single experiments where the rows are the genes and the columns are the samples and the rownames are the gene names. 
#' @param gene A gene name.
#' @param ylab The label for the y-axis.
#' @param legend_loc The position of the legend, default is "topright".
#' @param add_legend Logical. Should a legend of the names of the experiments be included.
#' @param add_main Logical. Should a header be included. Header would be the gene name.
#' @return A line plot.
#' @examples
#' plot_raw(gene = "SPBC2F12.09c")
#' plot_raw(data = profiles, gene = "SPBC2F12.09c", ylab = "profiles")

plot_raw <- function(data = single_exp, gene = "SPAC212.04c", 
                     ylab = "raw counts", legend_loc = "topright", 
                     add_legend = TRUE, add_main = FALSE) {
  
   main <- ifelse(add_main, gene, "")

   # extract gene from individual experiments
   gene_data <- c(data[[1]][gene, , drop = FALSE], 
                  data[[2]][gene, , drop = FALSE], 
                  data[[3]][gene, , drop = FALSE],
                  data[[4]][gene, , drop = FALSE], 
                  data[[5]][gene, , drop = FALSE], 
                  data[[6]][gene, , drop = FALSE])  
   yrange <- c(0, range(as.data.frame(gene_data), na.rm = TRUE)[2])
   
   matplot(c(0, 15, 30, 60, 120, 180), 
           t(as.data.frame(data[[1]][gene, , drop = FALSE])), 
           ty = "o", pch=19, ylab = ylab, xlab = "time [minutes]", 
           main = main, col = cols[1], ylim = yrange, cex.lab = 1.2)
   axis(3, at = c(0, 15, 30, 60, 120, 180), 
        labels = c("T0", "T1", "T2", "T3", "T4", "T5"))   
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[2]][gene, , drop = FALSE])), 
         col = cols[1], lty = 2, type = "o", pch = 19)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[3]][gene, , drop = FALSE])), 
         col = cols[1], lty = 3, type = "o", pch = 19)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[4]][gene, , drop = FALSE])), 
         col = cols[2], lty = 1, type = "o", pch = 19)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[5]][gene, , drop = FALSE])), 
         col = cols[2], lty = 2, type = "o", pch = 19)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[6]][gene, , drop = FALSE])), 
         col = cols[2], lty = 3, type = "o", pch = 19)
   
   if (add_legend)
     legend(legend_loc, col = cols[rep(1:2, each = 3)], lwd = 2, 
            lty = c(1, 2, 3, 1, 2, 3), 
            legend = c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3"))
}

########################################################
#' Plot mean profiles of single experiments
#'
#' @param data A list of mean profiles where the rows are the genes and the columns are the samples and the rownames are the gene names. 
#' @param gene A gene name.
#' @param legend_loc The position of the legend, default is "topright".
#' @param add_legend Logical. Should a legend of the names of the experiments be included.
#' @param add_main Logical. Should a header be included. Header would be the gene name.
#' @return A line plot.
#' @examples
#' plot_mprof(gene = "SPBC2F12.09c")

plot_mprof <- function(data = mprofiles, gene = "SPBC2F12.09c", 
                       legend_loc = "topright", add_legend = TRUE, 
                       add_main = FALSE) {
  
  main <- ifelse(add_main, gene, "")
  
  gene_data <- c(data[[1]][gene, , drop = FALSE], 
                 data[[2]][gene, , drop = FALSE])
                
  yrange <- c(0, range(as.data.frame(gene_data))[2])
  
  matplot(c(0, 15, 30, 60, 120, 180), 
          t(as.data.frame(data[[1]][gene, , drop = FALSE])), ty = "o",  
          ylab = "mean profiles", cex.lab = 1.2, xlab = "time [minutes]", 
          main = main, col = cols[1], pch = 19, ylim = yrange)
  axis(3, at = c(0, 15, 30, 60, 120, 180), 
       labels = c("T0", "T1", "T2", "T3", "T4", "T5"))   
  
  lines(c(0, 15, 30, 60, 120, 180), 
        t(as.data.frame(data[[2]][gene, , drop = FALSE])), 
        col = cols[2], lty = 1, type = "o", pch = 19)
  
  if (add_legend)
    legend(legend_loc, col = cols[1:2], lwd = 2, legend = c("WT", "Mut"))
}

########################################################
#' Plot transformed profiles of single experiments
#'
#' @param data A list of transformed profiles where the rows are the genes 
#' and the columns are the samples and the rownames are the gene names. 
#' @param gene A gene name.
#' @param ylab Label of the y-axis, typically the name of the transformation, 
#' e.g. "ALR".
#' @param legend_loc The position of the legend, default is "topright".
#' @param add_legend Logical. Should a legend of the names of the experiments 
#' be included.
#' @param add_main Logical. Should a header be included. Header would be the 
#' gene name.
#' @return A line plot.
#' @examples
#' plot_trafo(gene = "SPBC2F12.09c")

plot_trafo <- function(data = alr_object, gene = "SPBC2F12.09c", 
                       ylab = "ALR", legend_loc = "topright", 
                       add_legend = TRUE, add_main = FALSE) {
  
  main <- ifelse(add_main, gene, "")
  
  # get transformed data
  gene_data <- c(data[[1]][gene, , drop = FALSE], 
                 data[[2]][gene, , drop = FALSE])
  
  yrange <- c(0, range(as.data.frame(gene_data))[2])
  
  matplot(1:5, t(as.data.frame(data[[1]][gene, , drop = FALSE])), 
          ty = "o", pch = 19, ylab = ylab, cex.lab = 1.2, xaxt = "none", 
          xlab = "coordinates", main = main, col = cols[1], ylim = yrange)
  
  lines(1:5, t(as.data.frame(data[[2]][gene, , drop = FALSE])), 
        col = cols[2], lty = 1, type = "o", pch = 19)

  axis(1, at = 1:5, labels = c("C1", "C2", "C3", "C4", "C5"))
  
  if (add_legend)
    legend(legend_loc, col = cols[1:2], lwd = 2, legend = c("WT", "Mut"))
}

###############################################################
# visualization of the preprocessing steps to generate Figure 2

pdf("Figure2.pdf", width = 10, height = 6)
par(mfrow = c(2, 2), mar = c(4.1, 4.1, 2.1, 0.1), oma = c(0, 0, 3, 0))
plot_raw(data = single_exp, gene = "SPNCRNA.1164", add_main = FALSE)
plot_raw(data = profiles, gene = "SPNCRNA.1164", 
         add_main = FALSE, add_legend = TRUE, ylab = "profiles")
plot_mprof(add_main = FALSE, add_legend = TRUE, gene = "SPNCRNA.1164")
plot_trafo(add_main = FALSE, add_legend = TRUE, gene = "SPNCRNA.1164",
           ylab = "ALR transformed profiles")
mtext("SPNCRNA.1164", outer = TRUE, line = 1, cex = 2)
dev.off()

