library(HMMcopy)
rfile <- system.file("extdata", "normal.wig", package = "HMMcopy")
gfile <- system.file("extdata", "gc.wig", package = "HMMcopy")
mfile <- system.file("extdata", "map.wig", package = "HMMcopy")
normal_reads <- wigsToRangedData(rfile, gfile, mfile)
normal_reads[1000:1010, ]
normal_copy <- correctReadcount(normal_reads)
normal_copy[1000:1010, ]
copy = tumour_copy$copy
param <- data.frame(strength = 1e+07, e = 0.99, 
                    mu = quantile(tumour_copy$copy, na.rm = TRUE, prob = c(0.1, 
                                                                          0.25, 0.5, 0.75, 0.99)), lambda = 20, nu = 2.1, 
                    kappa = c(0.05,0.05,0.7, 0.1,  0.05) * 1000, 
                    m = 0, eta = c(5, 5,50, 5, 5) * 10000, gamma = 3, 
                    S = 0)
param$m <- param$mu
param$S <- ((sd(2^tumour_copy$copy[autosomes], na.rm = TRUE)/sqrt(nrow(param)))^2)
rownames(param) <- seq(1, 5)

chr <- tumour_copy$chr
chr <- as.factor(chr)
autosomes <- (chr != "X" & chr != "Y" & chr != "23" & 
      chr != "24" & chr != "chrX" & chr != "chrY" & chr != 
      "M" & chr != "MT" & chr != "chrM")
	  
maxiter=50
verbose = TRUE

K <- dim(param)[1]
N <- length(copy)
rho <- matrix(0, K, N)
py <- matrix(0, K, N)
mus <- matrix(0, K, maxiter)
lambdas <- matrix(0, K, maxiter)
kappa <- param$kappa
converged <- FALSE
Z <- rep(0, N)
Zcounts <- matrix(0, K, K)
loglik <- rep(0, maxiter)
chrs <- levels(chr)
chrsI <- vector("list", length(chrs))
for (i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
}
if (verbose) {
    message("Initialization")
}
i <- 1
mus[, 1] <- param$mu
lambdas[, 1] <- param$lambda
for (k in 1:K) {
    py[k, ] <- HMMcopy:::tdistPDF(copy, mus[k, i], lambdas[k, i], 
                        param$nu[k])
}
A <- matrix(0, K, K)
for (j in 1:K) {
    A[j, ] <- (1 - param$e[1])/(K - 1)
    A[j, j] <- param$e[1]
}
A_prior <- A
dirPrior <- A * param$strength[1]
loglik[i] <- -Inf


while (!converged && (i < maxiter)) {
    if (verbose) {
        message(paste("EM iteration:", i, "Log likelihood:", 
                      loglik[i]))
        message(converged)
    }
    i <- i + 1
    if (verbose) {
        message("Expectation")
    }
    Zcounts <- matrix(0, K, K)
    for (j in 1:length(chrsI)) {
        I <- chrsI[[j]]
		dyn.load("forward_backward.dll")
        output <- .Call("forward_backward", kappa, A, py[,I])
        rho[, I] <- output$rho
        loglik[i] <- loglik[i] + output$loglik
        Zcounts <- Zcounts + t(colSums(aperm(output$xi, 
                                             c(3, 2, 1))))
    }
    if (verbose) {
        message("Maximization")
    }
    mu_i <- mus[, i - 1]
    lambda_i <- lambdas[, i - 1]
    output <- HMMcopy:::estimateTNoiseParamsMap(copy[autosomes], mu_i, 
                                      lambda_i, param$nu, rho[, autosomes], param$eta, 
                                      param$m, param$gamma, param$S, param$kappa)
    mus[, i] <- output$mu_N
    lambdas[, i] <- output$lambda_N
    kappa <- output$pi
    for (k in 1:K) {
        py[k, ] <- HMMcopy:::tdistPDF(copy, mus[k, i], lambdas[k, 
                                                     i], param$nu[k])
    }
    priorA <- 0
    for (k in 1:K) {
        A[k, ] <- Zcounts[k, ] + dirPrior[k, ]
        A[k, ] <- HMMcopy:::normalize(A[k, ])
        priorA <- priorA + log(HMMcopy:::dirichletpdf(A_prior[k, ], 
                                            A[k, ]))
    }
    priorMu <- c()
    for (k in 1:K) {
        priorMu[k] <- log(dnorm(mus[k, i], param$mu[k], 
                                1))
    }
    loglik[i] <- loglik[i] + priorA + sum(priorMu)
    if (abs(loglik[i] - loglik[i - 1]) < 0.1 || loglik[i] < 
        loglik[i - 1]) {
        converged = 1
    }
}

if (converged) {
    i = i - 1
    if (verbose) {
        message("Re-calculating latest responsibilties for output")
    }
    for (j in 1:length(chrsI)) {
        I <- chrsI[[j]]
        output <- .Call("forward_backward", kappa, A, py[,I])
        rho[, I] <- output$rho
    }
}

if (verbose) {
    message("Optimal parameters found, segmenting and classifying")
}
segs <- vector("list", length(chrs))
for (c in 1:length(chrsI)) {
    I <- chrsI[[c]]
	dyn.load("viterbi.dll")
    output <- .Call("viterbi", log(kappa), log(A), log(py[,I]))
    Z[I] <- output$path
    segs[[c]] <- output$seg
}


