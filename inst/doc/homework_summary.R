## -----------------------------------------------------------------------------
x <- c(rep(1, 9),rep(2, 5),rep(3, 7),rep(4, 2),rep(5, 4))
b <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
a <- c("A", "B", "C", "D", "E")
d <- terrain.colors(5)
hist(x, breaks = b, labels = a, col = d)

## ----include=FALSE------------------------------------------------------------
x <- 1:1000; y <- log10(x)
b <- lm(y ~ x)
df <- summary(b)$coef
print(df)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(df)


## ---- fig.height=5------------------------------------------------------------
n <- 1e5
u <- runif(n)
sigma <- c(2, 5, 8, 11)
x1 <- sqrt(-2 * sigma[1] ^ 2 * log(1 - u)) 
x2 <- sqrt(-2 * sigma[2] ^ 2 * log(1 - u)) 
x3 <- sqrt(-2 * sigma[3] ^ 2 * log(1 - u)) 
x4 <- sqrt(-2 * sigma[4] ^ 2 * log(1 - u)) 
par(mfrow = c(2,2))
hist(x1, prob = TRUE, main = "sigma=2")
hist(x2, prob = TRUE, main = "sigma=5")
hist(x3, prob = TRUE, main = "sigma=8")
hist(x4, prob = TRUE, main = "sigma=11")

## ---- fig.height=6,eval=FALSE-------------------------------------------------
#  n <- 1000
#  x1 <- rnorm(n,0,1)
#  x2 <- rnorm(n,3,1)
#  r <- sample(c(0,1), n, replace = TRUE, prob = c(.75, .25))
#  r1 <- sample(c(0,1),n, replace = TRUE, prob = c(.4, .6))
#  r2 <- sample(c(0,1),n, replace = TRUE, prob = c(.6, .4))
#  r3 <- sample(c(0,1),n, replace = TRUE, prob = c(.5, .5))
#  y <- r * x1 + (1 - r) * x2
#  y1 <- r1 * x1 + (1 - r1) * x2
#  y2 <- r2 * x1 + (1 - r2) * x2
#  y3 <- r3 * x1 + (1 - r3) * x2
#  par(mfrow = c(2,2))
#  hist(y, col = 'light blue', prob = TRUE, main = "p1=0.75")
#  lines(density(y), col = "red", lwd = 3)
#  hist(y1, col = 'light blue', prob = TRUE, main = "p1=0.4")
#  lines(density(y1), col = "red", lwd = 3)
#  hist(y2, col = 'light blue', prob = TRUE, main = "p1=0.6")
#  lines(density(y2), col = "red", lwd = 3)
#  hist(y3, col = 'light blue', prob = TRUE, main = "p1=0.5")
#  lines(density(y3), col = "red", lwd = 3)
#  

## ---- fig.height=5------------------------------------------------------------
n <- 1e5
t <- 10
alpha <- 3; beta <- 4
for ( lambda in c(2, 4, 6, 8)) {
  Nt <- rpois(n, lambda * t)
  Xt <- rgamma(n, alpha * Nt, beta)
  E <- mean(Xt);  D <- var(Xt)
  E0 <- lambda * t * (alpha / beta)
  D0 <- lambda * t * (alpha * (alpha + 1)) / (beta ^ 2)
  EXt <- c(E, E0)
  VarX <- c(D, D0)
  v <- c("估计值", "理论值")
  a <- data.frame(v , EXt, VarX)
  b <- knitr::kable(a, digits = 3) 
  print(b)
}


## -----------------------------------------------------------------------------
MC.B <- function(x, m=10000){
  u <- runif(m)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
  g <- u ^ 2 * (1 - x[i] * u) ^ 2
  cdf[i] <- 30 * x[i] ^ 3 * mean(g)
  }
  cdf
}

x <- seq(.1, .9, .1)
mce <- MC.B(x)
pbe <- pbeta(x, 3, 3)
print(round(rbind(x,mce, pbe), 4))

## -----------------------------------------------------------------------------
MC.Ray <- function(x, R = 10000, sigma = 2, antithetic = TRUE){
  u <- runif(R/2)
  if(!antithetic) v <- runif(R/2) else
    v <- 1 - u
  u <- c(u, v)
  cdf <- numeric(length(x))
  for (i in 1:length(x)) {
    g <- x[i]^2 * u * exp(- (x[i] * u)^2 / (2 * sigma^2))
    cdf[i] <- mean(g) / sigma^2
  }
  cdf
}

m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.5
for (i in 1:m) {
  MC1[i] <- MC.Ray(x, R=1000, antithetic = FALSE)
  MC2[i] <- MC.Ray(x, R=1000)
}
print(c(sd(MC1), sd(MC2), (var(MC1) - var(MC2)) / var(MC1)))

## -----------------------------------------------------------------------------
m <- 10000
mc.e <- d <- numeric(2)
g <- function(x){
  exp(- x^2 / 2 + 2 * log(x)) /sqrt(2 * pi) * (x > 1)
}

#using f1
t <- rexp(m)
x <- t + 1
fg <- g(x) / exp(-t)
mc.e[1] <- mean(fg)
d[1] <- var(fg)

#using f2
t <- rnorm(m)
x <- abs(t) + 1
f2 <- sqrt(2 / pi) * exp(- t ^ 2 / 2)
fg <- g(x) / f2
mc.e[2] <- mean(fg)
d[2] <- var(fg)

print(round(rbind(mc.e, d), 5))

## -----------------------------------------------------------------------------
n <- 20
alpha <- .975
set.seed(123)
m <- 1000
LCL <- UCL <- numeric(m)
# 例6.4正态分布均值估计的置信区间覆盖率
for (i in 1:m) {
  x <- rnorm(n, mean = 0, sd = 2)
  UCL[i] <- mean(x)+sd(x)*qt(alpha, df= n-1)/sqrt(n)
  LCL[i] <- mean(x)-sd(x)*qt(alpha, df= n-1)/sqrt(n)
}
a <- mean(LCL<=0 & UCL>=0)
# 卡方分布均值估计的置信区间覆盖率
for (i in 1:m) {
  x <- rchisq(n, df = 2)
  UCL[i] <- mean(x)+sd(x)*qt(alpha, df= n-1)/sqrt(n)
  LCL[i] <- mean(x)-sd(x)*qt(alpha, df= n-1)/sqrt(n)
}
b <- mean(LCL<=2 & UCL>=2)
print(c(a,b))

## ----eval=FALSE---------------------------------------------------------------
#  n <- 20
#  alpha <- 0.05
#  set.seed(123)
#  m <- 1000
#  p1 <- p2 <- p3 <- numeric(m)
#  # (i)卡方分布
#  for (i in 1:m) {
#   x <- rchisq(n, df=1)
#   ttest <- t.test(x, alternative = "two.sided", mu = 1)
#   p1[i] <- ttest$p.value
#  }
#  p1.hat <- mean(p1<=alpha)
#  se1.hat <- sqrt(p1.hat *(1-p1.hat)/m)
#  # (ii)均匀分布
#  for (i in 1:m) {
#   x <- runif(n, min = 0, max = 2)
#   ttest <- t.test(x, alternative = "two.sided", mu = 1)
#   p2[i] <- ttest$p.value
#  }
#  p2.hat <- mean(p2<=alpha)
#  se2.hat <- sqrt(p2.hat *(1-p2.hat)/m)
#  # (iii)指数分布
#  for (i in 1:m) {
#   x <- rexp(n, 1)
#   ttest <- t.test(x, alternative = "two.sided", mu = 1)
#   p3[i] <- ttest$p.value
#  }
#  p3.hat <- mean(p3<=alpha)
#  se3.hat <- sqrt(p3.hat *(1-p3.hat)/m)
#  # 模拟结果（第一类错误率和估计的标准误差）
#  a <- c("(i)", "(ii)", "(iii)")
#  p <- c(p1.hat, p2.hat, p3.hat)
#  se <- c(se1.hat, se2.hat, se3.hat)
#  d <- data.frame(a, p, se)
#  c <- knitr::kable(d)
#  print(c)

## -----------------------------------------------------------------------------
n <- 500
alpha <- 0.05
set.seed(123)
m <- 1000
p1 <- p2 <- p3 <- numeric(m)
# (i)卡方分布
for (i in 1:m) {
 x <- rchisq(n, df=1) 
 ttest <- t.test(x, alternative = "two.sided", mu = 1)
 p1[i] <- ttest$p.value
}
p1.hat <- mean(p1<=alpha)
se1.hat <- sqrt(p1.hat *(1-p1.hat)/m)
# (ii)均匀分布
for (i in 1:m) {
 x <- runif(n, min = 0, max = 2) 
 ttest <- t.test(x, alternative = "two.sided", mu = 1)
 p2[i] <- ttest$p.value
}
p2.hat <- mean(p2<=alpha)
se2.hat <- sqrt(p2.hat *(1-p2.hat)/m)
# (iii)指数分布
for (i in 1:m) {
 x <- rexp(n, 1) 
 ttest <- t.test(x, alternative = "two.sided", mu = 1)
 p3[i] <- ttest$p.value
}
p3.hat <- mean(p3<=alpha)
se3.hat <- sqrt(p3.hat *(1-p3.hat)/m)
# 模拟结果（第一类错误率和估计的标准误差）
a <- c("(i)", "(ii)", "(iii)")
p <- c(p1.hat, p2.hat, p3.hat)
se <- c(se1.hat, se2.hat, se3.hat)
d <- data.frame(a, p, se)
c <- knitr::kable(d)
print(c)

## ---- eval=FALSE--------------------------------------------------------------
#  library(MASS)
#  #fist write a function to compute the sample Mardia's multivariate skewness statistic and test.
#  Mardia_sk <- function(X){
#    n <- nrow(X)
#    d <- ncol(X)
#    C <- X
#    for (j in 1:d) {
#      C[,j] <- X[,j] - mean(X[,j])
#    }
#    sigma_bar <- t(C)%*%C/n#compute the maximum likehihood estimator of covariance
#    A <-C %*%solve(sigma_bar)%*%t(C)
#    b <- sum(rowSums(A ^{3})) / (n^2)
#    T <- n*b/6
#    cv <- qchisq(.95, df=d*(d+1)*(d+2)/6)#critical value
#    as.integer(T > cv)
#  }
#  
#  set.seed(1234)
#  n <- c(10, 20, 30, 50, 100, 500)
#  mu <- c(0,0,0)
#  sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow = 3,ncol = 3)
#  p.reject <- numeric(length(n))
#  m <- 1000
#  for (i in 1:length(n)) {
#    p.reject[i] <- mean(replicate(m, expr = {
#      X <- mvrnorm(n[i], mu, sigma)
#      Mardia_sk(X)
#    }))
#  }
#  n <- as.character(n)
#  table <- rbind(n, p.reject)
#  knitr::kable(table ,digits = 4 ,align = "c")

## ---- fig.height=5, eval=FALSE------------------------------------------------
#  library(MASS)
#  Mardia_sk <- function(X){
#    n <- nrow(X)
#    d <- ncol(X)
#    C <- X
#    for (j in 1:d) {
#      C[,j] <- X[,j] - mean(X[,j])
#    }
#    sigma_bar <- t(C)%*%C/n#compute the maximum likehihood estimator of covariance
#    A <-C %*%solve(sigma_bar)%*%t(C)
#    b <- sum(rowSums(A ^{3})) / (n^2)
#    T <- n*b/6
#    cv <- qchisq(.95, df=d*(d+1)*(d+2)/6)#critical value
#    as.integer(T > cv)
#  }
#  
#  
#  set.seed(1234)
#  n <- 30
#  m <- 2500
#  epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
#  N <- length(epsilon)
#  pwr <- numeric(N)
#  mu1 <- mu2 <- c(0,0,0)
#  sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1), nrow = 3, ncol = 3)
#  sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100), nrow = 3, ncol = 3)
#  for (k in 1:N) {
#    e <- epsilon[k]
#    sktests <- numeric(m)
#    for (j in 1:m) {
#      index <- sample(c(1,2), size = n, replace = TRUE, prob = c(1-e, e))
#      X <- matrix(0, nrow = n, ncol = 3)
#      for (i in 1:n) {
#        if(index[i] == 1)
#          X[i,] <- mvrnorm(1, mu1, sigma1)
#        else
#          X[i,]<- mvrnorm(1, mu2, sigma2)
#      }
#      sktests[j] <-   Mardia_sk(X)
#    }
#    pwr[k] <- mean(sktests)
#  }
#  
#  plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylab = c(0,1))
#  abline(h = .05, lty = 3)
#  se <- sqrt(pwr*(1-pwr)/m)
#  lines(epsilon, pwr + se, lty=3)
#  lines(epsilon, pwr - se, lty=3)

## ---- eval=FALSE--------------------------------------------------------------
#  library(boot)
#  library(bootstrap)
#  set.seed(12345)
#  theta <- function(x,i){
#    lambda <- eigen(cov(x[i,]))$values
#    lambda[1] / sum(lambda)
#  }
#  obj <- boot(data = scor, statistic = theta, R = 2000)
#  original <- obj$t0
#  bias.boot <- mean(obj$t)-obj$t0
#  se.boot <- sd(obj$t)
#  table <- t(c(original, bias.boot, se.boot))
#  knitr::kable(table,align = "c", col.names = c("original", "bias.boot", "se.boot"))

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot)
set.seed(12345)
theta <- function(x,i){
  lambda <- eigen(cov(x[i,]))$values
  lambda[1] / sum(lambda)
}
n <- nrow(scor)
theta.hat <- theta(scor, 1:n)
theta.jack <- numeric(n)
for (i in 1:n) {
  theta.jack[i] <- theta(scor, (1:n)[-i])
}
bias.jack <- (n-1) * (mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1) * mean((theta.jack-theta.hat)^2))
table <- t(c(theta.hat, bias.jack, se.jack))
knitr::kable(table, align = "c", col.names = c("original", "bias.jack", "se.jack"))

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot)
set.seed(12345)
theta <- function(x,i){
  lambda <- eigen(cov(x[i,]))$values
  lambda[1] / sum(lambda)
}
obj <- boot(data = scor, statistic = theta, R = 2000)
ci <- boot.ci(obj, conf = 0.95, type = c("perc", "bca"))
print(ci)


## ---- eval=FALSE--------------------------------------------------------------
#  library(boot)
#  library(moments)
#  set.seed(12345)
#  boot.sk <- function(x,i){
#    skewness(x[i])
#  }
#  m <- 1000
#  n <- 20
#  ci.norm <- ci.basic <- ci.perc <- matrix(NA, m, 2)
#  cr.norm <- cr.basic <- cr.perc <- numeric(2)#the coverage rates
#  #normal population
#  for (i in 1:m) {
#    x <- rnorm(n,5,10)
#    a <- boot(data = x, statistic = boot.sk, R=1000)
#    ci <- boot.ci(a, conf = 0.95, type = c("norm", "basic","perc"))
#    ci.norm[i,] <- ci$norm[2:3]
#    ci.basic[i,] <- ci$basic[4:5]
#    ci.perc[i,] <- ci$percent[4:5]
#  }
#  cr.norm[1] <- mean(ci.norm[,1]<=0 & ci.norm[,2]>=0)
#  cr.basic[1] <- mean(ci.basic[,1]<=0 & ci.basic[,2]>=0)
#  cr.perc[1] <- mean(ci.perc[,1]<=0 & ci.perc[,2]>=0)
#  
#  #Chi-square distribution with 5 degrees of freedom
#  for (i in 1:m) {
#    x <- rchisq(n, df = 5)
#    sk0 <- sqrt(8/5)
#    a <- boot(data = x, statistic = boot.sk, R=1000)
#    ci <- boot.ci(a, conf = 0.95, type = c("norm", "basic","perc"))
#    ci.norm[i,] <- ci$norm[2:3]
#    ci.basic[i,] <- ci$basic[4:5]
#    ci.perc[i,] <- ci$percent[4:5]
#  }
#  cr.norm[2] <- mean(ci.norm[,1]<=sk0 & ci.norm[,2]>=sk0)
#  cr.basic[2] <- mean(ci.basic[,1]<=sk0 & ci.basic[,2]>=sk0)
#  cr.perc[2] <- mean(ci.perc[,1]<=sk0 & ci.perc[,2]>=sk0)
#  table <- data.frame(cr.norm, cr.basic, cr.perc,row.names = c("normal","chisq"))
#  rmarkdown::paged_table(table)

## -----------------------------------------------------------------------------
set.seed(122)
x <- rexp(15); y <- rexp(15)
cat("x:",x,'\n');cat("y:",y,'\n')
R <- 999
z <- c(x,y)
K <- 1:30
reps <- numeric(R)
t0 <- cor.test(x,y, method = "spearman")$estimate
for (i in 1:R) {
  k <- sample(K, size =15 , replace = FALSE)
  xstar <- z[k]
  ystar <- z[-k]
  reps[i] <- cor(xstar,ystar, method = "spearman")
}
p <- mean(abs(c(t0, reps)) >= abs(t0))
round(c(p, cor.test(x,y)$p.value), 4)

## ---- eval=FALSE--------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(Ball)
#  library(energy)
#  library(MASS)
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z,0);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1)
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#    (i1 + i2) / (k * n)
#  }
#  
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R, sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)
#  }

## ---- eval=FALSE--------------------------------------------------------------
#  mu1 <- c(0,0,0)
#  sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
#  mu2 <- c(0.5,-0.5,0.5)
#  sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
#  n1=n2=20
#  n <- n1+n2
#  N = c(n1,n2)
#  k=3
#  R=999
#  m=100
#  set.seed(1234)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    mydata1 <- mvrnorm(n1,mu1,sigma1)
#    mydata2 <- mvrnorm(n2,mu2,sigma2)
#    mydata <- rbind(mydata1,mydata2)
#    p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ---- eval=FALSE--------------------------------------------------------------
#  mu1 <- c(0,0,0)
#  sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
#  mu2 <- c(0.5,-0.5,0.5)
#  sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
#  n1=n2=20
#  n <- n1+n2
#  N = c(n1,n2)
#  k=3
#  R=999
#  m=100
#  set.seed(1234)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    mydata1 <- mvrnorm(n1,mu1,sigma1)
#    mydata2 <- mvrnorm(n2,mu2,sigma2)
#    mydata <- rbind(mydata1,mydata2)
#    p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ---- eval=FALSE--------------------------------------------------------------
#  n1=n2=20
#  n <- n1+n2
#  N = c(n1,n2)
#  k=3
#  R=999
#  m=100
#  set.seed(1234)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    mydata1 <- as.matrix(rt(n1,1,2),ncol=1)
#    mydata2 <- as.matrix(rt(n2,2,5),ncol=1)
#    mydata <- rbind(mydata1,mydata2)
#    p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ---- eval=FALSE--------------------------------------------------------------
#  n1=n2=20
#  n <- n1+n2
#  N = c(n1,n2)
#  k=3
#  R=999
#  m=100
#  set.seed(1234)
#  p.values <- matrix(NA,m,3)
#  rbimodel<-function(n,mu1,mu2,sd1,sd2){
#    index=sample(1:2,n,replace=TRUE)
#    x=numeric(n)
#    index1<-which(index==1)
#    x[index1]<-rnorm(length(index1), mu1, sd1)
#    index2<-which(index==2)
#    x[index2]<-rnorm(length(index2), mu2, sd2)
#    return(x)
#  }
#  for(i in 1:m){
#    mydata1 <- as.matrix(rbimodel(n1,0,0,1,2),ncol=1)
#    mydata2 <- as.matrix(rbimodel(n2,1,1,4,3),ncol=1)
#    mydata <- rbind(mydata1,mydata2)
#    p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ---- eval=FALSE--------------------------------------------------------------
#  mu1 <- c(0,0,0)
#  sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
#  mu2 <- c(0.5,-0.5,0.5)
#  sigma2 <- matrix(c(2,0,0,0,2,0,0,0,2),nrow=3,ncol=3)
#  n1=10
#  n2=100
#  n <- n1+n2
#  N = c(n1,n2)
#  k=3
#  R=999
#  m=100
#  set.seed(1234)
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    mydata1 <- mvrnorm(n1,mu1,sigma1)
#    mydata2 <- mvrnorm(n2,mu2,sigma2)
#    mydata <- rbind(mydata1,mydata2)
#    p.values[i,1] <- eqdist.nn(mydata,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(mydata,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=mydata1,y=mydata2,num.permutations=R,seed=i*2846)$p.value
#  }
#  alpha <- 0.05;
#  pow <- colMeans(p.values<alpha)
#  pow

## ---- fig.height=4------------------------------------------------------------
set.seed(12345)
fcau<- function(x,theta,eta){
  stopifnot(theta > 0)
  return( 1 / (theta * pi *(1 + ((x - eta) / theta)^2)))
}
theta <- 1; 
eta <- 0
n <- 1e4
x <- numeric(n)
sigma <- 2
x[1] <-rnorm(1, mean=0, sd=sigma)
k <- 0
u <- runif(n)
for (i in 2:n) {
  xt <- x[i-1]
  y <- rnorm(1, mean=xt, sd=sigma)
  sup <- fcau(y, theta, eta) * dnorm(xt, mean=y, sd=sigma)
  sub <- fcau(xt, theta, eta) * dnorm(y, mean=xt, sd=sigma)
  if( u[i] <= sup/sub) {
    x[i] <- y
  } else{
    x[i] <- xt
    k <- k+1
  }
}

cat(" the rate of rejected candidate points=",k/n)
par(mfrow=c(1,2))
b <- 1001   #discard the burn-in sample
y <- x[b:n]
p <- ppoints(10)
QR <- qcauchy(p)  #quantiles of Cauchy
Q <- quantile(x, p)
qqplot(QR, Q, main="", xlab="Cauchy Quantiles", ylab="Sample Quantiles")
abline(0,1,col='blue',lwd=2)
hist(y, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR, fcau(QR, theta, eta))

## -----------------------------------------------------------------------------
set.seed(666)
M <- 5000
burn <- 1000
X <- matrix(0,nrow = M, ncol = 2)
n <- 10
a <- 1
b <- 2
X[1,] <- c(0, .5)
for (i in 2:M) {
  x2 <- X[i-1, 2]
  X[i,1] <- rbinom(1, n, x2)
  x1 <- X[i,1]
  alpha <- x1+a
  beta <- n-x1+b
  X[i,2] <- rbeta(1, alpha, beta)
}

b <- burn+1
x <- X[b:M, ]
cm <- colMeans(x)
result1 <- c("colmean", cm[1], cm[2])
knitr::kable(t(result1), align = "c", col.names = c("", "X1","X2"))
result2 <- data.frame(cov(x),row.names = c("X1","X2"))
rmarkdown::paged_table(result2)
result3 <- data.frame(cor(x),row.names = c("X1","X2") )
rmarkdown::paged_table(result3)
plot(x, main = "", cex = 0.5, xlab = bquote(X[1]), ylab = bquote(X[2]), ylim = range(x[,2]))

## ---- fig.height=4------------------------------------------------------------
G.R.M <- function(psi){
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, "var")
  w <- mean(psi.w)
  v.hat <- w*(n-1)/n + (B/n)
  r.hat <- v.hat / w
  return(r.hat)
}

fcau<- function(x){
  return(1 / ( pi *(1 + x^2)))
}

cauchy.chain <- function(sigma, N, X1){
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rnorm(1, xt, sigma)
    sup <- fcau(y) * dnorm(xt, y, sigma)
    sub <- fcau(xt) * dnorm(y, xt, sigma)
    r <- sup/sub
    if(u[i] <= r) x[i] <- y else
      x[i] <- xt
  }
  return(x)
}
set.seed(123)
sigma <- sqrt(2)
k <- 4
n <- 15000
b <- 500
x0 <- c(-10,-5,5,10)
X <- matrix(0, nrow = k, ncol = n)
for (i in 1:k) {
  X[i,] <- cauchy.chain(sigma, n, x0[i])
}
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) {
  psi[i,] <- psi[i,]/(1:ncol(psi))
}
for (i in 1:k) {
  if(i==1){
    plot((b+1):n,psi[i, (b+1):n],ylim=c(-2,2), type="l", xlab='Index', ylab=bquote(phi))
    }
  else{
    lines(psi[i, (b+1):n], col=i)
    }
}
rhat <- rep(0, n)
for (j in (b+1):n){
  rhat[j] <- G.R.M(psi[,1:j])
}
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
G.R.M <- function(psi){
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)
  B <- n * var(psi.means)
  psi.w <- apply(psi, 1, "var")
  w <- mean(psi.w)
  v.hat <- w*(n-1)/n + (B/n)
  r.hat <- v.hat / w
  return(r.hat)
}
set.seed(356)
a <- 1;b <- 2;n <- 10
f.chain <- function(x0,y0,M){
  Z <- matrix(0,M,2)
  Z[1,1] <- x0
  Z[1,2] <- y0
  for (i in 2:M) {
    Z[i,1] <- rbinom(1, n, Z[i-1,2])
    alpha <- Z[i,1]+a
    beta <- n-Z[i,2]+b
    Z[i,2] <- rbeta(1,alpha,beta)
  }
  return(Z)
}
m <- 15000
k <- 4
x0 <- c(1,2,3,4)
y0 <- c(.2, .3, .4, .5)
b <- 1000
X <- matrix(0, nrow = k, ncol = m)
Y <- matrix(0, nrow = k, ncol = m)
for (i in 1:k) {
  z <- f.chain(x0[i],y0[i],m)
  X[i,] <- z[,1]
  Y[i,] <- z[,2]
}
psix <- t(apply(X, 1, cumsum))
psiy <- t(apply(Y, 1, cumsum))
# for X
for (i in 1:nrow(psix)) {
  psix[i,] <- psix[i,]/(1:ncol(psix))
}
for (i in 1:k) {
  if(i==1){
    plot((b+1):m,psix[i, (b+1):m],ylim=c(-.02,.02), type="l", xlab='Index', ylab=bquote(phi))
    }
  else{
    lines(psix[i, (b+1):m], col=i)
    }
}
rhatx <- rep(0, m)
for (j in (b+1):m){
  rhatx[j] <- G.R.M(psix[,1:j])
}
plot(rhatx[(b+1):m], type="l", xlab="X", ylab="R")
abline(h=1.2, lty=2)
# for Y
for (i in 1:nrow(psiy)) {
  psiy[i,] <- psiy[i,]/(1:ncol(psiy))
}
for (i in 1:k) {
  if(i==1){
    plot((b+1):m,psiy[i, (b+1):m],ylim=c(-0.005,.005), type="l", xlab='Index', ylab=bquote(phi))
    }
  else{
    lines(psiy[i, (b+1):m], col=i)
    }
}
rhaty <- rep(0, m)
for (j in (b+1):m){
  rhaty[j] <- G.R.M(psiy[,1:j])
}
plot(rhaty[(b+1):m], type="l", xlab="Y", ylab="R")
abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
# write a function to compute the kth term
kth <- function(a,k){
  d <- length(a)
  ad <- numeric(d)
  for (j in 1:d) {
    ad[j] <- a[j]^2
  }
  euc2 <- sum(ad)
  b <- euc2^(k+1) / ((2*k+1)*(2*k+2))
  m <- exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+1+d/2))
  return((-1)^k * b * m / (factorial(k)* (2^k)))
}
#modify the function so that it computes and returns the sum
fl <- function(a,n){
  result <- numeric(n)
  for (k in 1:n) {
    result[k] <- kth(a, k-1)
  }
  return(sum(result))
}
#Evaluate the sum when a=(1,2)
a <- c(1,2)
n <- c(4,10,21,51,101)
result <- numeric(length(n))
for (i in 1:length(n)) {
  result[i] <- fl(a, n[i])
}
N <- as.character(n)
table <- data.frame(N,result)
knitr::kable(t(table), align = "c")

## ---- eval=FALSE--------------------------------------------------------------
#  # 11.4 find the intersection points A(k)
#  sk1 <- function(a,k){
#    ck <- sqrt(a^2 * k / (k+1- a^2))
#    pt(ck, df=k)
#  }
#  f1 <- function(a,k){sk1(a,k)-sk1(a,k-1)}
#  ak1 <- function(k){
#    solve <- uniroot(function(a){f1(a,k)}, lower = 0.5, upper = 2)
#    return(solve$root)
#  }
#  
#  # 11.5 Write a function to solve the equation for a
#  sk2 <- function(a,k){
#    ck <- sqrt(a^2 * k / (k+1- a^2))
#    res <- pt(ck, df=k)-pt(-ck, df=k)
#    return(res)
#  }
#  f2 <- function(a,k){sk2(a,k)-sk2(a,k-1)}
#  ak2 <- function(k){
#    solve <- uniroot(function(a){f2(a,k)}, lower = 0.5, upper = 2)
#    return(solve$root)
#  }
#  k <- c(4:25, 100, 500, 1000)
#  root <- matrix(0, nrow = 3, ncol = length(k))
#  root[1, ] <- k
#  for (i in 1:length(k)) {
#    root[2,i] <- ak1(k[i])
#    root[3,i] <- ak2(k[i])
#  }
#  
#  row.names(root) <- c("k","A(k)_11.4","A(k)_11.5")
#  knitr::kable(t(root), align = "c")

## ---- eval=FALSE--------------------------------------------------------------
#  # Use the E-M algorithm to estimate lambda
#  y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
#  a <- sum(y)
#  n <- length(y)
#  m=0
#  for (i in 1:n) {
#    if(y[i] == 1)
#      m <- m+1
#    else
#      m <- m
#  }
#  #the result with the observed data MLE
#  lambda_bar <- a /(n-m)
#  # Use the E-M algorithm to estimate lambda
#  lambda <- 2 # initial estimated value for lambda
#  N <- 2000
#  lambda.old <- lambda + 1
#  tol <- .Machine$double.eps ^ 0.5
#  for (k in 1:N) {
#    lambda <- (a + m*lambda) / n
#    if(abs(lambda-lambda.old) / lambda.old < tol) break
#    lambda.old <- lambda
#  }
#  print(list(lambda.EM=lambda,lambda.MLE=lambda_bar ,iter=k, tol=tol))
#  

## ---- eval=FALSE--------------------------------------------------------------
#  trims <- c(0, 0.1, 0.2, 0.5)
#  x <- rcauchy(100)
#  a <- lapply(trims, function(trim) mean(x, trim = trim))
#  b <- lapply(trims, mean, x = x)
#  knitr::kable(data.frame(unlist(a),unlist(b)), align = "c")

## ---- eval=FALSE--------------------------------------------------------------
#  attach(mtcars)
#  formulas <- list(
#  mpg ~ disp,
#  mpg ~ I(1 / disp),
#  mpg ~ disp + wt,
#  mpg ~ I(1 / disp) + wt
#  )
#  # with a for loop
#  myout3 <- vector("list", length(formulas))
#  for (i in seq_along(formulas)) {
#    myout3[[i]] <- lm(formulas[[i]], mtcars)
#  }
#  # with lapply
#  lq3 <- lapply(formulas, function(f) lm(formula = f, data = mtcars))
#  rsq <- function(mod) summary(mod)$r.squared
#  #储存数据并制作表格
#  tt <- matrix(0, nrow = length(formulas), ncol = 2)
#  tt[,1] <- unlist(lapply(myout3, rsq))
#  tt[,2] <- unlist(lapply(lq3, rsq))
#  colnames(tt) <- c("for loop", "lapply")
#  knitr::kable(t(tt), align = "c", col.names = c("mpq~disp","mpg~I(1/disp)","mpq~disp+wt","mpg~I(1/disp)+wt"))

## ---- eval=FALSE--------------------------------------------------------------
#  bootstraps <- lapply(1:10, function(i) {
#  rows <- sample(1:nrow(mtcars), rep = TRUE)
#  mtcars[rows, ]
#  })
#  # with a for loop
#  myout4 <- vector("list", length(bootstraps))
#  for (i in seq_along(bootstraps)) {
#    myout4[[i]] <- lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp)
#  }
#  # with lapply
#  lq4 <- lapply(bootstraps, function(x) lm(formula = mpg~disp, data = x))
#  rsq <- function(mod) summary(mod)$r.squared
#  #储存数据并制作表格
#  tt <- matrix(0, nrow = length(bootstraps), ncol = 2)
#  tt[,1] <- unlist(lapply(myout4, rsq))
#  tt[,2] <- unlist(lapply(lq4, rsq))
#  colnames(tt) <- c("for loop", "lapply")
#  knitr::kable(t(tt), align = "c", col.names = as.character(1:10))
#  

## ---- eval=FALSE--------------------------------------------------------------
#  #a) Compute the standard deviation of every column in a numeric data frame.
#  lq <- data.frame(x1=runif(20,0,2), x2=rnorm(20,0,1), x3=rchisq(20,8), x4=rexp(20, 2))
#  a <- vapply(lq, sd, FUN.VALUE =numeric(1))
#  knitr::kable(t(round(a,4)), align = "c")
#  
#  #b) Compute the standard deviation of every numeric column in a mixed data frame.
#  name <- c("Alice","Lvy","May","Emma","Edith")
#  gender <- c("女","女","男","女","男")
#  score1 <- runif(5, 70,100)
#  score2 <- runif(5, 80,90)
#  score3 <- runif(5, 60,100)
#  b <- data.frame(x1=name, x2=gender, x3=score1, x4=score2, x5=score3)
#  # 构造函数筛选数值型的列数据并计算标准差，非数值型返回NA
#  sd.num <- function(i,b){
#    c <- vapply(b, class, character(1))
#    if(c[i] == "numeric")
#      sd(b[,i])
#    else return(NA)
#  }
#  res <- vapply(seq_along(b),function(i) sd.num(i=i, b=b), numeric(1))
#  knitr::kable(t(res), align = "c", col.names = names(b))

## ---- eval=FALSE--------------------------------------------------------------
#  # 构造mcsapply函数
#  mcsapply <- function(X,FUN) {
#    library(parallel)
#    cl.cores <- detectCores()
#    cl <- makeCluster(cl.cores)
#    g<-  parSapply(cl,X, FUN)
#    stopCluster(cl)
#    return(g)
#  }
#  
#  #实验一下是否加快运算速度
#  test <- replicate(1e5, t.test(rchisq(30, df=1)), simplify = FALSE)
#  a <- system.time(mcsapply(test,function(x){unlist(x)[3]}))
#  b <- system.time(sapply(test,function(x) {unlist(x)[3]}))
#  list(a,b)

## ---- warning=FALSE, eval=FALSE-----------------------------------------------
#  library(Rcpp)
#  # Write an Rcpp function for Exercise 9.8
#  sourceCpp(code='
#  #include <cmath>
#  #include <Rcpp.h>
#  using namespace Rcpp;
#  
#  //[[Rcpp::export]]
#    NumericMatrix gibbsC(int n, int a, int b, int M){
#      NumericMatrix mat(M, 2);
#      double x1=0;
#      double x2=0.5;
#      double alpha=0;
#      double beta=0;
#      for(int i=0 ; i < M; i++){
#        x1 = rbinom(1, n, x2)[0];
#        alpha = x1 + a;
#        beta = n - x1 + b;
#        x2 = rbeta(1, alpha, beta)[0];
#        mat(i,0) = x1;
#        mat(i,1) = x2;
#      }
#      return(mat);
#    }
#  ')
#  
#  # Write an R function for Exercise 9.8
#  gibbsR <- function(n,a,b,M){
#    X <- matrix(0,nrow = M, ncol = 2)
#    X[1,] <- c(0, .5)
#    for (i in 2:M) {
#      x2 <- X[i-1, 2]
#      X[i,1] <- rbinom(1, n, x2)
#      x1 <- X[i,1]
#      alpha <- x1+a
#      beta <- n-x1+b
#      X[i,2] <- rbeta(1, alpha, beta)
#    }
#    return(X)
#  }
#  
#  #Compare the corresponding generated random numbers
#  set.seed(9903)
#  n <- 20;a <- 2; b <- 4
#  M <- 3500
#  burn <- 500
#  c <- burn + 1
#  randomnum.R <- gibbsR(n,a,b,M)[c:M, ]
#  randomnum.C <- gibbsC(n,a,b,M)[c:M, ]
#  par(mfrow=c(1,2))
#  qqplot(randomnum.R[,1], randomnum.C[,1],main="", xlab = "X with pure R language", ylab = "x with Rcpp function")
#  abline(a=0,b=1,col='blue', lwd=2)
#  qqplot(randomnum.R[,2], randomnum.C[,2],main="", xlab = "Y with pure R language", ylab = "Y with Rcpp function")
#  abline(a=0,b=1,col='blue', lwd=2)
#  
#  #Compare the computation time
#  library(microbenchmark)
#  times <- microbenchmark(tm.R=gibbsR(n,a,b,M),tm.c=gibbsC(n,a,b,M))
#  table <- summary(times)
#  knitr::kable(table, align = "c")

