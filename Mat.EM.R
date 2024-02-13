#' Check validity of object for initialization
#'
#' @param initial A list of a fitted mixture of matrix-variate
#'     Gaussians.
#' @param dim A numeric vector of length three indicating the data
#'     dimensions.
#' @param tolerance: A numeric >= 0. Differences smaller than
#'     'tolerance' for the row sums of 'gamma' ar not reported.
#' @return A logical; TRUE if all checks succeed.

check_initial <- function(initial, dim, tol = sqrt(.Machine$double.eps)) {
    stopifnot(is.list(initial))
    stopifnot(is.vector(dim), is.numeric(dim), length(dim) == 3)
    stopifnot(is.vector(tol), is.numeric(tol), length(tol) == 1,
              tol > 0)
    stopifnot(all(c("gamma", "Mu", "Psi", "Sigma") %in% names(initial)))
    
    stopifnot(is.matrix(initial$gamma), is.numeric(initial$gamma),
              nrow(initial$gamma) == dim[3])
    stopifnot(all(initial$gamma >= 0),
              all(abs(rowSums(initial$gamma) - 1) < tol))
    
    G <- ncol(initial$gamma)
    stopifnot(is.array(initial$Mu), is.numeric(initial$Mu),
              all(dim(initial$Mu) == c(dim[1:2], G)))
    stopifnot(is.array(initial$Psi), is.numeric(initial$Psi),
              all(dim(initial$Psi) == c(rep(dim[2], 2), G)))
    stopifnot(is.array(initial$Sigma), is.numeric(initial$Sigma),
              all(dim(initial$Sigma) == c(rep(dim[1], 2), G)))
    TRUE
}

#' Density of the matrix-variate normal distribution
#'
#' @param X A numeric array of matrix-variate observations of dimension p x d x n.
#' @param Mu A numeric matrix of dimension p x d.
#' @param Psi A numeric matrix of dimension d x d.
#' @param Sigma A numeric matrix of dimension p x p.
#' @param log A logical; if TRUE, probabilities p are given as log(p).
#' @return A vector of length 'n' giving the density values. 
dmatnorm <- function(X, Mu, Psi, Sigma, log = FALSE) {
    stopifnot(is.array(X), is.numeric(X), length(dim(X)) == 3)
    p <- dim(X)[1]
    d <- dim(X)[2]
    stopifnot(is.matrix(Mu), is.numeric(Mu), dim(Mu) == c(p, d))
    stopifnot(is.matrix(Psi), is.numeric(Psi), dim(Psi) == c(d, d))    
    stopifnot(is.matrix(Sigma), is.numeric(Sigma), dim(Sigma) == c(p, p))
    stopifnot(is.vector(log), is.logical(log), length(log) == 1)
    
    Z <- sweep(X, 1:2, Mu, "-")
    Psic <- tryCatch(base::chol(Psi), error = function(e) e)
    Sigmac <- tryCatch(base::chol(Sigma), error = function(e) e)
    if (inherits(Psic, "error") || inherits(Sigmac, "error")) {
        X.is.mu <- apply(abs(Z), 3, sum) == 0
        logretval <- rep.int(-Inf, dim(X)[3])
        logretval[X.is.mu] <- Inf
    }
    else {
        rss <- apply(Z, 3, function(z) {
            tmp <- backsolve(Sigmac, z, transpose = TRUE)
            sum(diag(tcrossprod(backsolve(Psic, t(tmp), transpose = TRUE))))
        })
        logretval <- - d * sum(log(diag(Sigmac))) - p * sum(log(diag(Psic))) -
            0.5 * d * p * log(2 * pi) - 0.5 * rss
    }
    if (log) 
        logretval
    else exp(logretval)
}

## Internal function from package flexmix.
log_row_sums <- function(m) {
    M <- m[cbind(seq_len(nrow(m)), max.col(m, "first"))]
    M + log(rowSums(exp(m - M)))
}

## Internal function from package flexmix.
printIter <- function(iter, logLik, label = "Log-likelihood") 
cat(formatC(iter, width = 4), label, ":", formatC(logLik, width = 12, 
    format = "f"), "\n")

#' EM algorithm to fit a mixture of matrix-variate Gaussians
#'
#' @param Y A numeric array of matrix-variate observations of
#'     dimension p x d x n.
#' @param initial A fitted mixture for initialization, e.g., obtained
#'     using 'MatTrans.EM'.
#' @param row_model Character indicating the structure of the row-wise
#'     covariance.
#' @param col_model Character indicating the structure of the
#'     column-wise covariance.
#' @param control List specifying 'tolerance',
#'     'iter.max" and 'verbose' for the EM algorithm.
#' @return A list similar to a single fitted model returned by
#'     'MatTrans.EM' as an element in 'result'.
Mat.EM <- function(Y, initial, row_model = c("diag", "full"),
                   col_model = c("diag", "full", "ar1"), control = NULL) {
    stopifnot(is.array(Y), length(dim(Y)) == 3)
    check_initial(initial, dim(Y))
    
    row_model <- match.arg(row_model)
    col_model <- match.arg(col_model)
    default_control <- list(tolerance = 1e-5, iter.max = 100, verbose = TRUE)
    if (is.null(control)) {
        control <- default_control
    }  
    stopifnot(is(control, "list"), length(control) > 0)
    stopifnot(length(control) == length(names(control)))
    stopifnot(all(names(control) %in% names(default_control)))
    default_control[names(control)] <- control
    control <- default_control
    
    postscaled <- initial$gamma
    Mu <- initial$Mu
    Psi <- initial$Psi
    Sigma <- initial$Sigma

    llh <- -Inf
    converged <- FALSE
    G <- ncol(postscaled)
    n <- nrow(postscaled)
    p <- nrow(Sigma)
    d <- nrow(Psi)

    Psic <- array(0, dim(Psi))
    Sigmac <- array(0, dim(Sigma))

    pars <- G * d * p +
        G * switch(row_model,
                   diag = p,
                   full = p * (p + 1) / 2) +
        G * switch(col_model,
                   diag = d,
                   full = d * (d + 1) / 2,
                   ar1 = 2 * d - 1) - 1
    for (iter in 1:control$iter.max) {
        ## M-step
        tau <- colMeans(postscaled)
        n_g <- colSums(postscaled)
        for (g in 1:G) {
            Mu[,, g] <- 
                apply(sweep(Y, 3, postscaled[, g], "*"), 1:2, sum) / n_g[g]
            Z <- sweep(Y, 1:2, Mu[,,g])
            Psic <- base::chol(Psi[,,g])
            Sigma[,,g] <- Reduce("+", lapply(1:n, function(i)
                postscaled[i, g] * crossprod(backsolve(Psic, t(Z[,,i]), transpose = TRUE)))) / (d * n_g[g])
            if (row_model == "diag") {
                Sigma[,,g] <- diag(diag(Sigma[,,g]))
            }
            Sigmac <- base::chol(Sigma[,,g])
            Psi[,,g] <- Reduce("+", lapply(1:n, function(i)
                postscaled[i, g] * crossprod(backsolve(Sigmac, Z[,,i], transpose = TRUE)))) / (p * n_g[g])
            if (col_model == "diag") {
                Psi[,,g] <- diag(diag(Psi[,,g]))
            } else if (col_model == "ar1") {
                T_g <- diag(nrow = d)
                T_g[cbind(2:d, 1:(d-1))] <- - Psi[,,g][cbind(2:d, 1:(d-1))] / diag(Psi[,,g])[-d]
                D_g <- T_g %*% Psi[,,g] %*% t(T_g)
                T_g_inv <- forwardsolve(T_g, diag(d))
                Psi[,,g] <- tcrossprod(sweep(T_g_inv, 2, sqrt(diag(D_g)), "*"))
            }
            alpha <- det(Psi[,,g])^(1/d)
            Psi[,,g] <- Psi[,,g] / alpha
            Sigma[,,g] <- Sigma[,,g] * alpha
        }
        ## E-step
        postunscaled <- do.call("cbind", lapply(1:G, function(g)
            dmatnorm(Y, Mu[,,g], Psi[,,g], Sigma[,,g], log = TRUE)))
        postunscaled <- sweep(postunscaled, 2, log(tau), "+")
        logpostunscaled <- postunscaled
        postunscaled <- exp(postunscaled)
        postscaled <- exp(logpostunscaled - log_row_sums(logpostunscaled))
        llh.old <- llh
        llh <- sum(log_row_sums(logpostunscaled))
        if (abs(llh - llh.old)/(abs(llh) + 0.1) < control$tolerance) {
            if (control$verbose > 0) {
                printIter(iter, llh)
                cat("converged\n")
            }
            converged <- TRUE
            break
        }
        if (control$verbose && (iter %% control$verbose == 0)) 
            printIter(iter, llh)
    }
    return(list(tau = tau, Sigma = Sigma, Psi = Psi, Mu = Mu, gamma = postscaled,
                iter = iter, pars = pars, id = max.col(postscaled), flag = converged,
                ll = llh, bic = -2 * llh + pars * log(n)))
}

