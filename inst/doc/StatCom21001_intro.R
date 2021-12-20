## ---- eval=FALSE--------------------------------------------------------------
#  MDS <- function(X, n){
#    m <- nrow(X)
#    B <- matrix(0, nrow = m, ncol = m)
#    D <- matrix(0, nrow = m, ncol = m)
#    for (i in 1:m) {
#      for (j in 1:m) {
#        D[i,j] <- sum((X[i,]-X[j,])^2)
#      }
#    }
#    Di <- colSums(D) / m
#    Dj <- rowSums(D) / m
#    Dij <- sum(D) / (m^2)
#    for (i in 1:m) {
#      for (j in 1:m) {
#        B[i,j] <- (D[i,j] + Dij - Di[i]-Dj[j]) / (-2)
#      }
#    }
#    B_eigen_values <- eigen(B)$values
#    B_eigen_vectors <- eigen(B)$vectors
#    a <- sort(B_eigen_values, decreasing = TRUE)[1:n]
#    index <- numeric(n)
#    for (k in seq_along(a)) {
#      index[k] <- which(B_eigen_values==a[k], arr.ind = TRUE)
#    }
#    v <- diag(B_eigen_values[index])
#    u <- B_eigen_vectors[, index]
#    z <- sqrt(v) %*% t(u)
#    return(z)
#  }

## ---- eval=FALSE--------------------------------------------------------------
#   X <- matrix(c(3,5,2,3,4,8,9,6,1,7,5,9,8,6,4,4),nrow=4)
#   MDS(X, n=2)

## ---- eval=FALSE--------------------------------------------------------------
#  para_EM <- function(x, mean, sd=NULL){
#    num <- length(mean)
#    if(is.null(sd)){
#      sd <- rep(1, num)
#      }
#    epsilon <- 1e-4
#    probs <- rep(1/num, num)
#    mu_s <- mean
#    sigma_s <- sd ^ 2
#    n <- length(x)
#    while (TRUE) {
#      es <- matrix(0, nrow = n, ncol = num)
#      for (j in seq(num)) {
#        es[,j] <- probs[j] * dnorm(x, mean=mu_s[j], sd=sqrt(sigma_s[j]))
#      }
#      es <- es / rowSums(es)
#      sigma_s_p <- sigma_s
#      for (j in seq(num)) {
#        sigma_s[j] <- sum(es[,j] * (x-mu_s[j]) ^ 2) / sum(es[,j])
#        mu_s[j] <- sum(x * es[,j]) / sum(es[,j])
#        probs[j] <- mean(es[,j])
#      }
#      if(max(abs(sigma_s_p-sigma_s)) < epsilon){ break}
#  
#    }
#    a <- list(mu=mu_s, sd=sqrt(sigma_s), prob=probs)
#    return(a)
#  }
#   x <- rnorm(1000,1,1)
#   para_EM(x, mean = c(0,1), sd= c(1,1))

## ---- eval=FALSE--------------------------------------------------------------
#  n <- 1000
#  mean_s <- c(1, 7)
#  y <- sample(c("head", "tail"), size = n, replace = TRUE, prob = c(0.25, 0.75))
#  x <- rnorm(n = 1000, mean = mean_s[1])
#  tails <- y %in% c("tail")
#  x[tails] <- rnorm(sum(tails), mean = mean_s[2])
#  para_EM(x, mean = c(0, 1), sd = c(1, 1))

