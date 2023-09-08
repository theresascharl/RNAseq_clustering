require("tidyverse")
require("robCompositions")
require("flexclust")
require("RColorBrewer")
require("ggpubr")


## get posteriors and calculate dbsi ----
get.si <- function(M = M, genes = NULL, type = "3way") {
  
  type <- match.arg(type, choices = c("3way", "2way"))
  
  ## get posteriors 
  if(type == "2way") pos <- M$z
  else   pos <- M$best.result[[1]]$gamma
  
  ## number of classes
  n <- ncol(pos)
  
  ## get largest and second largest posteriors incl classes
  res <- matrix(nrow=nrow(pos), ncol=8)
  rownames(res) <- genes
  colnames(res) <- c("p1", "p2", "pac", "c1", "ca", "r", "dbsi", "co")
  res[,1] <- apply(pos, 1, max)
  res[res[,1]>0.99999,1] <- 0.99999
  res[,2] <- apply(pos, 1, function(x) sort(x, partial = n-1)[n-1])
  res[res[,2]<0.00001,2] <- 0.00001
  res[,3] <- 1 - 2*(res[,2]/(res[,1]+res[,2]))
  res[,4] <- apply(pos, 1, which.max)
  res[,5] <- apply(pos, 1, function(x) order(x, decreasing = TRUE)[2])
  res[,6] <- log(res[,1]/res[,2])
  res[,7] <- res[,6]/max(res[,6]) #dbsi
  res <- as.data.frame(res)
  
  mean_returned <- mean(res[,7])
  
  # average dbsi per cluster
  ave <- tapply(res$dbsi, res$c1, mean)
  
  # order clusters by decreasing dbsi
  o <- order(ave, decreasing = TRUE)
  
  for(i in 1:n) res[res$c1 == o[i],8] <- i
  
  # now order within clusters by dbsi
  ores <- res[order(res$co),]
  for(i in 1:n)
    ores[ores$co==i,] <- ores[ores$co==i,][order(ores$dbsi[ores$co==i],
                                                 decreasing = TRUE),] 
  
  return(list(result = ores, mean = mean_returned, order = o))
}

## silhouette plot after Raymaekers in R package classmap ----
silplot_pos <- function (result, classCols = NULL, showLegend = TRUE, 
                         showClassNumbers = FALSE, showCases = FALSE, 
                         drawLineAtAverage = FALSE, method = "dbsi",
                         topdown = TRUE, main = NULL, summary = TRUE) 
{
  result$co <- factor(result$co)
  nlab <- length(unique(result$co))
  
  lvls <- 1:nlab
  if (is.null(classCols)) {
    classCols <- hcl.colors(nlab)
  }
  else {
    if (!is.vector(classCols)) {
      stop("\n classCols should be a vector")
    }
    if (length(classCols) < nlab) {
      stop(paste0("\n The number of classCols should be at", 
                  " least length(vcrout$levels) = ", nlab, "."))
    }
    classCols <- classCols[seq_len(nlab)]
  }
  
  yintv <- result$co
  indsv <- which(!is.na(yintv))
  ayint <- sort(unique(yintv))
  
  #if(method=="PAC") si <- result$pac
  #else
  si <- result$dbsi
  
  if (is.null(main)) {
    main <- ""
  }
  df <- as.data.frame(cbind(class = yintv, si, name = yintv))
  if (topdown) {
    df <- df[order(-df$class, df$si), ]
  }
  else {
    df <- df[order(df$class, -df$si), ]
  }
  avswidth <- mean(df$si)
  df$name <- factor(df$name, levels = levels(ayint))
  df$name2 <- 1:nrow(df)
  df$name2 <- as.factor(df$name2)
  df$class <- as.factor(df$class)
  if (showClassNumbers) {
    mapping <- aes_string(x = "name2", y = "si", color = "class", 
                          fill = "class")
  }
  else {
    df$label <- factor(df$name, levels = levels(ayint))
    mapping <- aes_string(x = "name2", y = "si", color = "label", 
                          fill = "label")
  }
  gg <- ggplot(df, mapping) + 
    geom_bar(stat = "identity", show.legend = showLegend, size = 0.05, 
             width = 0.75) + 
    labs(y = paste0("dbsi"), 
         caption = paste0("\nOverall average dbsi: ", 
                          round(avswidth, 2)), x = "", title = paste0(main))
  
  silwidths_perclass <- round(sapply(lvls, function(y) 
    mean(df$si[which(df$label == y)])), 2)
  n <- tapply(df$class, df$class, length)

  templabels <- c()
  for (i in seq_len(length(unique(df$class)))) {
    sbar <- sprintf("%1.2f", silwidths_perclass[i])
    templabel <- paste("paste(bar(s), \" = \", ", "\"", 
                       sprintf("%1.2f", silwidths_perclass[i]), "\"", ", \" \", ", 
                       "\"", n[i], "\"", ")")
    templabels <- c(templabels, templabel)
  }
  gg <- gg + scale_fill_manual(values = classCols[ayint], 
                               aesthetics = c("colour", "fill"), labels = parse(text = templabels), 
                               name = "Cluster dbsi and size") + theme(legend.text.align = 0)
  if (topdown) 
    gg <- gg + coord_flip()
  if (drawLineAtAverage) {
    gg <- gg + geom_hline(yintercept = avswidth, linetype = "dashed", 
                          color = "red")
  }
  # if (!showCases) {
  #   if (topdown) {
  #     gg <- gg + theme(axis.text.y = element_blank(), 
  #                      axis.ticks.y = element_blank())
  #   }
  #   else {
  #     gg <- gg + theme(axis.text.x = element_blank(), 
  #                      axis.ticks.x = element_blank())
  #   }
  # }
  # else {
  #   if (topdown) {
  #     gg <- gg + theme(axis.text.y = element_text(angle = 0))
  #   }
  #   else {
  #     gg <- gg + theme(axis.text.x = element_text(angle = 90))
  #   }
  # }
  gg <- gg + theme(plot.title = element_text(hjust = 0.5))
  gg <- gg + theme(axis.line.x = element_line(size = 0.25))
  gg <- gg + theme(axis.ticks.y = element_blank())
  cln_pos <- round(tapply(as.numeric(df$name2), df$class, mean))
  gg <- gg + scale_x_discrete(expand = c(0, 0.013 * length(indsv)),
                              labels = as.character(1:nlab),
                              breaks = as.character(cln_pos))
  gg <- gg + scale_y_continuous(expand = c(0, 0), limits = c(0,1))
  gg <- gg + theme(panel.background = element_blank())
  gg <- gg + theme(plot.title.position = "plot")
  gg <- gg + theme(axis.title = element_text(size = 10))
  gg <- gg + theme(plot.caption = element_text(size = 10, 
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
cluster_map <- function(M = M, genes = genes, type = "3way", 
                        mprofiles = mprofiles, npcol = 4) 
{
  res <- get.si(M, genes = genes, type=type)
  
  type <- match.arg(type, choices = c("3way", "2way"))
  
  if(type == "3way")  {
    
    M1 <- M$best.result[[1]]
    
    ## estimated means reordered by average dbsi
    mu <- M1$Mu[,,res$order]
    
    ## means of Exp1 and Exp2 
    nexp <- dim(mu)[1]
    ntimes <- dim(mu)[2] + 1
    
    ## cluster centers in transformed space ----
    trafo_mu <- lapply(seq(dim(mu)[1]), function(x) t(mu[x, ,]))
    
    ## class assignment ordered by average dbsi ----
    cl <- res$result$co
  }
  
  else
  {
    nexp <- length(mprofiles)
    ntimes <- dim(mprofiles[[1]])[2]
    mu <- t(M$parameters$mean)[res$order,]
    trafo_mu <- list()
    for(i in 1:nexp)
      trafo_mu[[i]] <- mu[,(1:(ntimes-1))+(i-1)*(ntimes-1)]
   ## class assignment ----
    cl <- res$result$co
  }
  
  ## cluster centers in original space ----
  orig_mu <- lapply(trafo_mu, function(x) 
    addLRinv(x, cnames = paste("T",1:ntimes,sep="_"), ivar = 1))
  
  ## mean profiles in orig space using subset of genes in the order of the dbsi ----
  Exp <- lapply(mprofiles, function(x) x[rownames(res$result),])
  
  ## matrix of the same number of rows containing the coordinates of the centers ----
  centers <- lapply(orig_mu, function(x) x[cl,])
  
  ## combine by flattening out the experiments ----
  Exp_all <- matrix(unlist(Exp), ncol= ntimes*nexp)
  cent <- matrix(unlist(centers), ncol= ntimes*nexp)
  
  ## calculate the Euclidean distance between genes and centers ----
  d <- rep(NA, length(cl))
  for(i in 1:length(d)) d[i] <- as.numeric(distEuclidean(Exp_all[i,,drop=FALSE],
                                                         cent[i,,drop=FALSE]))
  dmax <- max(d)
  
  ## reorder the dbsi to order of genes as used for the clustering ----
  ## not needed anymore res1 <- res$result[genes,]
  res1 <- res$result
  
  ## divide the distance by the max distance
  mat1 <- data.frame(cluster = cl, dist = d/dmax, dbsi = res1$dbsi)
  rownames(mat1) <- rownames(res1)

  ## plot dbsi vs dist for complete data as a panel plot
  k <- max(cl)
  mat1$cluster <- factor(mat1$cluster, labels = paste("cluster", 1:k))
  
  require(ggplot2)
    gg <- ggplot(mat1,aes(x=dist,y=dbsi,color=cluster))+
    scale_color_manual(values=hcl.colors(k))+
    scale_x_continuous(name = "farness to cluster center", 
                       limits=c(-0.1, 1.09)) +
    scale_y_continuous(limits=c(-0.1, 1.09)) +
    geom_point()+coord_fixed()+
    stat_chull(fill = NA)+ theme_bw()+
    facet_wrap(~cluster, ncol=npcol)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position = "none")
    
  print(gg)
  return(invisible(mat1))
}  

## add chull in ggplot2 ----
## code taken from 
## https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
