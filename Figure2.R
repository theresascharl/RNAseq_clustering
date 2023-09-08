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
           ty = "l", ylab = ylab, xlab = "time [minutes]", 
           main = main, col = cols[1], ylim = yrange, cex.lab = 1.2)
   
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[2]][gene, , drop = FALSE])), 
         col = cols[1], lty = 2)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[3]][gene, , drop = FALSE])), 
         col = cols[1], lty = 3)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[4]][gene, , drop = FALSE])), 
         col = cols[2], lty = 1)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[5]][gene, , drop = FALSE])), 
         col = cols[2], lty = 2)
   lines(c(0, 15, 30, 60, 120, 180), 
         t(as.data.frame(data[[6]][gene, , drop = FALSE])), 
         col = cols[2], lty = 3)
   
   if (add_legend)
     legend(legend_loc, col = cols[rep(1:2, each = 3)], lwd = 2, 
            lty = c(1, 2, 3, 1, 2, 3), 
            legend = c("WT1", "WT2", "WT3", "Mut1", "Mut2", "Mut3"))
}

plot_raw(gene = "SPBC2F12.09c")

plot_raw(data = profiles, gene = "SPBC2F12.09c", ylab = "profiles")

########################################################
plot_mprof <- function(data = mprofiles, gene = "SPBC2F12.09c", 
                       legend_loc = "topright", add_legend = TRUE, 
                       add_main = FALSE) {
  
  main <- ifelse(add_main, gene, "")
  
  gene_data <- c(data[[1]][gene, , drop = FALSE], 
                 data[[2]][gene, , drop = FALSE])
                
  yrange <- c(0, range(as.data.frame(gene_data))[2])
  
  matplot(c(0, 15, 30, 60, 120, 180), 
          t(as.data.frame(data[[1]][gene, , drop = FALSE])), ty = "l",  
          ylab = "mean profiles", cex.lab = 1.2, xlab = "time [minutes]", 
          main = main, col = cols[1], ylim = yrange)
  
  lines(c(0, 15, 30, 60, 120, 180), 
        t(as.data.frame(data[[2]][gene, , drop = FALSE])), 
        col = cols[2], lty = 1)
  
  if (add_legend)
    legend(legend_loc, col = cols[1:2], lwd = 2, legend = c("WT", "Mut"))
}

plot_mprof()

########################################################
plot_trafo <- function(data = alr_object, gene = "SPBC2F12.09c", 
                       ylab = "ALR", legend_loc = "topright", 
                       add_legend = TRUE, add_main = FALSE) {
  
  main <- ifelse(add_main, gene, "")
  
  # get transformed data
  gene_data <- c(data[[1]][gene, , drop = FALSE], 
                 data[[2]][gene, , drop = FALSE])
  
  yrange <- c(0, range(as.data.frame(gene_data))[2])
  
  matplot(1:5, t(as.data.frame(data[[1]][gene, , drop = FALSE])), 
          ty = "l", ylab = ylab, cex.lab = 1.2, xaxt = "none", 
          xlab = "coordinates", main = main, col = cols[1], ylim = yrange)
  
  lines(1:5, t(as.data.frame(data[[2]][gene, , drop = FALSE])), 
        col = cols[2], lty = 1)

  axis(1, at = 1:5, labels = c("C1", "C2", "C3", "C4", "C5"))
  
  if (add_legend)
    legend(legend_loc, col = cols[1:2], lwd = 2, legend = c("WT", "Mut"))
}

plot_trafo()

###################################################
# visualization of the preprocessing steps 

pdf("Figure2.pdf", width = 10, height = 6)
par(mfrow = c(2, 2), mar = c(4.1, 4.1, 0.1, 0.1), oma = c(0, 0, 3, 0))
plot_raw(data = single_exp, gene = "SPNCRNA.1164", add_main = FALSE)
plot_raw(data = profiles, gene = "SPNCRNA.1164", 
         add_main = FALSE, add_legend = FALSE, ylab = "profiles")
plot_mprof(add_main = FALSE, add_legend = FALSE, gene = "SPNCRNA.1164")
plot_trafo(add_main = FALSE, add_legend = FALSE, gene = "SPNCRNA.1164",
           ylab = "ALR transformed profiles")
mtext("SPNCRNA.1164", outer = TRUE, line = 1, cex = 2)
dev.off()

