rm(list = ls()) ; cat("\014")

fr <- function(y, x, betapm) {      
  
  p <- length(betapm)/2
  b <- sapply(1:p, function(j) betapm[j] - betapm[p+j] )
  
  as.numeric(crossprod(y - x %*% b))
}

gr <- function(y, x, betapm){
  
  p <- length(betapm)/2
  xy <- crossprod(x, y)
  xx <- crossprod(x, x)
  bp <- betapm[1:p]
  bn <- betapm[(p+1):(2*p)]
  
  as.numeric(
  rbind(-2 * xy + 2 * xx %*% (bp - bn),
         2 * xy - 2 * xx %*% (bp - bn)))
  
}


set.seed(1)
p <- 10
x = mvtnorm::rmvnorm(1000, mean = rep(0, p), sigma = diag(p))
beta <-  rep(0, p)
beta[c(1,3)] <- 2
y <- as.numeric(x %*% beta + rnorm(1000))

t <- sum(abs(beta)) + 10
A <- rbind(rep(-1, 2*p),
           diag(2*p))
c <- c(-t, rep(0, 2*p))

cop <- constrOptim(theta = rep(0.00001, 2*p), f = fr, grad = gr, y = y, x = x, 
                  ui = A, ci = c)

hat <- cop$par
round(sapply(1:p, function(j) hat[j] - hat[p+j] ), 5)

glm <- glmnet::cv.glmnet(x = x, y = y)
coef(glm)
