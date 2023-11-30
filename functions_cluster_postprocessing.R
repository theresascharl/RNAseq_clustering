library("tidyverse")
library("robCompositions")
library("flexclust")
library("RColorBrewer")
library("ggpubr")


#' Get posteriors and calculate dbsi 
#'
#' @param M A cluster result either from package MatTransMix or Mclust. 
#' @param genes A vector of gene names corresponding to the rows in the 
#' cluster result.
#' @param type Either "2way" (Mclust) or "3way" (MatTransMix) cluster object.
#' @return A list containing the data.frame 'result' with 8 variables: 
#'    p1 largest posterior (where max is set to 0.99999)
#'    p2 second largest posterior (where min is set to 0.00001)
#'    pac probability of the alternative class after Raymaekers in 
#'        R package classmap
#'    c1 cluster number with largest posterior
#'    ca cluster number with second largest posterior
#'    r log ratio between p1 and p2
#'    dbsi ratio of r and max(r)
#'    co newly assigned cluster number where cluster 1 has to largest dbsi
#'    
#' the 'order' of the original clusters and
#' the 'mean' of the dbsi values as an indicator for the overall cluster quality.
#' @examples
#' dbsi_3way <- get_dbsi(M, genes = genes)

get_dbsi <- function(M = M, genes = NULL, type = "3way") {
  
  type <- match.arg(type, choices = c("3way", "2way"))
  
  ## get posteriors 
  pos <- if (type == "2way") M$z
         else M$best.result[[1]]$gamma
  
  ## number of classes
  n <- ncol(pos)
  
  ## get largest and second largest posteriors incl classes
  res <- matrix(nrow = nrow(pos), ncol = 8)
  rownames(res) <- genes
  colnames(res) <- c("p1", "p2", "pac", "c1", "ca", "r", "dbsi", "co")
  res[, 1] <- apply(pos, 1, max)
  res[, 1] <- pmin(res[, 1], 0.99999) 
  res[, 2] <- apply(pos, 1, function(x) sort(x, partial = n-1)[n-1])
  res[, 2] <- pmax(res[, 2], 0.00001)
  res[, 3] <- 1 - 2 * (res[, 2] / (res[, 1] + res[, 2]))
  res[, 4] <- apply(pos, 1, which.max)
  res[, 5] <- apply(pos, 1, function(x) order(x, decreasing = TRUE)[2])
  res[, 6] <- log(res[, 1] / res[, 2])
  res[, 7] <- res[, 6] / max(res[, 6]) # dbsi
  res <- as.data.frame(res)
  
  mean_returned <- mean(res[, 7])
  
  ## average dbsi per cluster
  ave <- tapply(res$dbsi, res$c1, mean)
  
  ## order clusters by decreasing dbsi
  o <- order(ave, decreasing = TRUE)
  
  for (i in 1:n) res[res$c1 == o[i], 8] <- i
  
  ## now order within clusters by dbsi
  ores <- res[order(res$co), ]
  for (i in 1:n) {
      ores[ores$co==i, ] <- ores[ores$co == i,][order(ores$dbsi[ores$co == i],
                                                      decreasing = TRUE), ]
  }
  
  return(list(result = ores, mean = mean_returned, order = o))
}

## dbsi plot after Raymaekers in R package classmap ----

#' Density-based silhouette information plot 
#'
#' @param result The 'result' element of the returned object of function
#' 'get_dbsi'. 
#' @param showLegend Logical. Should a legend of the cluster dbsi and size be
#' included.
#' @param drawLineAtAverage Logical. Should average dbsi lines be included.
#' @param topdown Logical. Should the genes be ordered topdown or bottomup.
#' @param main Logical. Should a header be included.
#' @param summary Logical. Should a summary be printed in the console.
#' @return A dbsi information plot.
#' @examples
#' dbsi_plot(dbsi_3way$result)

dbsi_plot <- function (result, showLegend = TRUE, drawLineAtAverage = FALSE, 
                       topdown = TRUE, main = NULL, summary = TRUE) 
{
  result$co <- factor(result$co)
  nlab <- nlevels(result$co)
  
  lvls <- 1:nlab
  classCols <- hcl.colors(nlab)

  yintv <- result$co
  indsv <- which(!is.na(yintv))
  ayint <- sort(unique(yintv))
  
  si <- result$dbsi
  
  if (is.null(main)) {
    main <- ""
  }
  df <- data.frame(class = as.integer(levels(yintv))[as.integer(yintv)],
                   si = si, name = yintv)
  if (topdown) {
    df <- df[order(-df$class, df$si), ]
  } else {
    df <- df[order(df$class, -df$si), ]
  }
  avswidth <- mean(df$si)
  df$name <- factor(df$name, levels = levels(ayint))
  df$name2 <- 1:nrow(df)
  df$name2 <- as.factor(df$name2)
  df$class <- as.factor(df$class)
  df$label <- factor(df$name, levels = levels(ayint))
  
  gg <- ggplot(df, aes(x = name2, y = si, colour = class, fill = class)) + 
    geom_bar(stat = "identity", show.legend = showLegend, 
             width = 0.75) + 
    labs(y = paste0("dbsi"), 
         caption = paste0("\nOverall average dbsi: ", 
                          round(avswidth, 2)), x = "", title = paste0(main))
  
  silwidths_perclass <- round(sapply(lvls, function(y) 
      mean(df$si[which(df$label == y)])), 2)
  n <- tapply(df$class, df$class, length)

  templabels <- vapply(seq_len(length(unique(df$class))), function(i) {
      sbar <- sprintf("%1.2f", silwidths_perclass[i])
      paste("paste(bar(s), \" = \", ", "\"", 
            sprintf("%1.2f", silwidths_perclass[i]), "\"", ", \" \", ", 
            "\"", n[i], "\"", ")")
  }, character(1))
  gg <- gg +
      scale_fill_manual(values = classCols[ayint], 
                        aesthetics = c("colour", "fill"),
                        labels = parse(text = templabels), 
                        name = "Cluster dbsi and size") +
      theme(legend.text.align = 0)
  if (topdown) 
      gg <- gg +
          coord_flip()
  if (drawLineAtAverage) {
      gg <- gg +
          geom_hline(yintercept = avswidth, linetype = "dashed", 
                     color = "red")
  }
  gg <- gg +
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.line.x = element_line(linewidth = 0.25)) + 
      theme(axis.ticks.y = element_blank())
  cln_pos <- round(tapply(as.numeric(df$name2), df$class, mean))
  gg <- gg +
      scale_x_discrete(expand = c(0, 0.013 * length(indsv)),
                       labels = as.character(1:nlab),
                       breaks = as.character(cln_pos)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
      theme(panel.background = element_blank()) + 
      theme(plot.title.position = "plot") + 
      theme(axis.title = element_text(size = 10)) + 
      theme(plot.caption = element_text(size = 10, 
                                        hjust = 0, 
                                        margin = ggplot2::margin(0, 0, 0, 0)))
  if (summary) {
      ave <- tapply(df$si, df$class, mean)
      n <- tapply(df$class, df$class, length)
      sil.sum <- data.frame(classNumber = names(ave), 
                            classLabel = lvls[as.numeric(names(ave))], 
                            classSize = n, classAveSi = round(ave, 2), 
                            stringsAsFactors = TRUE)
      print(sil.sum, row.names = FALSE)
  }
  return(gg)
}


## add compactness of clusters

#' Cluster map 
#'
#' @param M A cluster result either from package MatTransMix or Mclust. 
#' @param genes A vector of gene names corresponding to the rows in the 
#' cluster result.
#' @param type Either "2way" (Mclust) or "3way" (MatTransMix) cluster object.
#' @param mprofiles A list of mean profiles where the rows are the genes and the columns are the samples and the rownames are the gene names.
#' @param npcol Number of columns used in the panel plot.
#' @return A list of 2 containing the dataframe 'info' with
#'   rows corresponding to the clustered genes and variables
#'   cluster the hard cluster membership
#'   dist the Euclidean distance to the cluster center standardized by dividing
#'   by the max distance
#'   dbsi the density-based silhouette information and
#' A panel plot where for each cluster the dbsi as a measure of cluster 
#' separation in the transformed space is plotted against the distance of genes 
#' to their cluster center in the original space.
#' @examples
#' cluster_map(M = M, genes = genes, type = "3way", mprofiles = mprofiles, 
#' npcol = 4)$plot
#' 
cluster_map <- function(M = M, genes = genes, type = c("3way", "2way"),
                        mprofiles = mprofiles, npcol = 4) 
{
  res <- get_dbsi(M, genes = genes, type = type)
  type <- match.arg(type)
  
  if (type == "3way")  {
    M1 <- M$best.result[[1]]
    
    ## estimated means reordered by average dbsi
    mu <- M1$Mu[, , res$order]
    
    ## means of Exp1 and Exp2 
    nexp <- dim(mu)[1]
    ntimes <- dim(mu)[2] + 1
    
    ## cluster centers in transformed space ----
    trafo_mu <- lapply(seq(dim(mu)[1]), function(x) t(mu[x, ,]))
    
    ## class assignment ordered by average dbsi ----
    cl <- res$result$co
  } else {
    nexp <- length(mprofiles)
    ntimes <- dim(mprofiles[[1]])[2]
    mu <- t(M$parameters$mean)[res$order, ]
    trafo_mu <- list()
    for (i in 1:nexp) {
        trafo_mu[[i]] <- mu[, (1:(ntimes-1))+(i-1)*(ntimes-1)]
    }
    ## class assignment ----
    cl <- res$result$co
  }
  
  ## cluster centers in original space ----
  orig_mu <- lapply(trafo_mu, function(x) 
    addLRinv(x, cnames = paste("T", 1:ntimes, sep = "_"), ivar = 1))
  
  ## mean profiles in orig space using subset of genes in the order of the dbsi ----
  Exp <- lapply(mprofiles, function(x) x[rownames(res$result), ])
  
  ## matrix of the same number of rows containing the coordinates of the centers ----
  centers <- lapply(orig_mu, function(x) x[cl, ])
  
  ## combine by flattening out the experiments ----
  Exp_all <- matrix(unlist(Exp), ncol= ntimes*nexp)
  cent <- matrix(unlist(centers), ncol= ntimes*nexp)
  
  ## calculate the Euclidean distance between genes and centers ----
  d <- vapply(seq_along(cl), function(i)
      as.numeric(distEuclidean(Exp_all[i ,, drop = FALSE], cent[i, , drop = FALSE])),
      numeric(1))
  dmax <- max(d)
  
  ## reorder the dbsi to order of genes as used for the clustering ----
  ## not needed anymore res1 <- res$result[genes,]
  res1 <- res$result
  
  ## divide the distance by the max distance
  info <- data.frame(cluster = cl, dist = d/dmax, dbsi = res1$dbsi)
  rownames(info) <- rownames(res1)

  ## plot dbsi vs dist for complete data as a panel plot
  k <- max(cl)
  info$cluster <- factor(info$cluster, labels = paste("cluster", 1:k))
  
  require("ggplot2")
  gg <- ggplot(info, aes(x = dist, y = dbsi,color = cluster)) +
      scale_color_manual(values = hcl.colors(k)) +
      scale_x_continuous(name = "distance to cluster center", 
                         limits = c(-0.1, 1.09)) +
      scale_y_continuous(limits = c(-0.1, 1.09)) +
      geom_point() +
      coord_fixed() +
      stat_chull(fill = NA) + 
      theme_bw() +
      facet_wrap(~cluster, ncol = npcol) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(legend.position = "none")
  return(list(info = info, plot = gg))
}  

## add chull in ggplot2 ----
## code taken from 
## https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                         data[chull(data$x, data$y), , drop = FALSE]
                     },
                     required_aes = c("x", "y"))

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' Calculate the ICL from posteriors of a MatTransMix object 
#'
#' @param object A cluster result from package MatTransMix. 
#' @return The ICL of a cluster result. 
#' @examples
#' icl_M(M)

icl_M <- function(object, ...)
{
  z <- object$best.result[[1]]$gamma
  n <- nrow(z)
  if (is.null(z)) {
    z <- matrix(1, nrow = n, ncol = 1)
  }
  C <- matrix(0, n, ncol(z))
  for (i in 1:n) {
    C[i, which.max(z[i,])] <- 1
  }
  object$best.bic + 2 * sum(C * ifelse(z > 0, log(z), 0))
}
