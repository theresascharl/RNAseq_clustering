library("robCompositions")
library("flexclust")
library("flexmix")
library("RColorBrewer")
palette(brewer.pal(8, "Dark2"))

## generate artificial data
n <- 100

## a data set with 4 multivariate normal groups in Euclidean space ----
set.seed(123)
myNclus <- rbind(mvtnorm::rmvnorm(n, mean = c(-3, -3),
                                  sigma = diag(rep(0.1, 2))),
                 mvtnorm::rmvnorm(n, mean = c(4, 0),
                                  sigma = diag(1:2)),
                 mvtnorm::rmvnorm(1.5 *n, mean = c(-2, 4),
                                  sigma = diag(c(1.8, 1))),
                 mvtnorm::rmvnorm(2 *n, mean = c(1, 1),
                                  sigma = matrix(c(0.1, 0.09, 0.09, 0.1), 2)))
colnames(myNclus) <- c("x1", "x2")

## transformation to the simplex
data1 <- addLRinv(myNclus, cnames = c("x1", "x2", "x3"), ivar = 1)

## a data set with 5 multivariate normal groups in Euclidean space ----
set.seed(123)
myNclus2 <- rbind(mvtnorm::rmvnorm(2*n, mean = c(-5, -6)),
                  mvtnorm::rmvnorm(n, mean = c(0, 0),
                                   sigma = diag(rep(0.1, 2))),
                  mvtnorm::rmvnorm(2 *n, mean = c(-1, -1),
                                   sigma = matrix(c(1, -0.9, -0.9, 1), 2)),
                  mvtnorm::rmvnorm(n, mean = c(7, 5),
                                   sigma = diag(rep(0.1, 2))),
                  mvtnorm::rmvnorm(n, mean = c(0, 5),
                                   sigma = diag(rep(0.1, 2))))
colnames(myNclus2) <- c("x1", "x2")

## transformation to the simplex
data2 <- addLRinv(myNclus2, cnames = c("x1", "x2", "x3"), ivar = 1)

pdf("Figure1.pdf")
par(mfrow = c(2, 2), mar = c(2.1, 3.1, 2, 0.2))
col1 <- rep(brewer.pal(4, "Dark2"), c(100, 100, 150, 200))
col2 <- alpha(col1, alpha = 0.4)
ternaryDiag(data1, col = col2, pch = 19, main = "Simplex")
plot(myNclus, col = col2, xlab = "", ylab = "", xlim = c(-10, 10),
     ylim = c(-10, 10), pch = 19, main = "Euclidean Space")
col1 <- rep(brewer.pal(5, "Dark2"), c(200, 100, 200, 100, 100))
col2 <- alpha(col1, alpha = 0.4)
ternaryDiag(data2, col = col2, pch = 19, main = "Simplex")
plot(myNclus2, col = col2, xlab = "", ylab = "", xlim = c(-10, 10),
     ylim = c(-10, 10), pch = 19, main = "Euclidean Space")
dev.off()

