dgp.screening <- function(n, share.lc = 0.1, reporting.bias = FALSE){
  #' simulated data of a LC screening trial
  ### 1) Creating covariates ----
  x1 <- 50 + rgamma(n, shape = 3.5, scale = 2.8)                                      # age (50+)
  x2 <- rbinom(n, 1, 0.16)                                                            # sex (~16% females)
  x3 <- sample(1:7, n, replace = TRUE, prob = c(0.1, 0.2, 0.15, 0.15, 0.1, 0.2, 0.1)) # education
  x4 <- sample(1:5, n, replace = TRUE, prob = c(0.05, 0.1, 0.6, 0.2, 0.05))           # health
  x5 <- rnorm(n, 27, 3)                                                               # bmi
  x6 <- 10 * (1 + rweibull(n, shape = 1.8, scale = 1)) ; x6[x6<0] <- 0                # num. cigarettes
  
  ## Smoking years:
  start.age                  <- rnorm(n, 16, 1)
  ex.smokers                 <- sample(1:n, floor(n*0.5)) # Let 50% of the participants be ex-smokers
  draw                       <- rnorm(length(ex.smokers), 20, 1)
  yrs.since.stop             <- rep(0, n)
  yrs.since.stop[ex.smokers] <- draw
  x7                         <- x1 - start.age - yrs.since.stop                       # smk years
  x7[x7 < 0] <- 0
  
  #Intermediate lifeyears:
  #set 5 lifeyears as base survival
  Base_lifeyears = 5
  #Change lifeyears based on x covariates: higher age reduces LY, female sex increases, education increases, health increases, BMI below 20 and above 25 decrease
  #number of cigarettes decrease, current smoking decreases, time smoked increases, time since smoking cessation decreases.
  Lifeyear_coeffs = c(-0.1,2,0.05,0.1,-0.2,-0.01,0.2,0.01,-0.02)
  Lifeyears_intermediate    <-  Base_lifeyears +  (Lifeyear_coeffs[1]*x1 ) + (Lifeyear_coeffs[2]*x2 ) + (Lifeyear_coeffs[3]*x3 ) + (Lifeyear_coeffs[4]*x4 ) + (Lifeyear_coeffs[5]*ifelse(x5<20 |x5>25,1,0 )) + (Lifeyear_coeffs[6]*x6 ) + (Lifeyear_coeffs[7]*x7 )
  Lifeyears_intermediate <-ifelse(Lifeyears_intermediate<=0,0,Lifeyears_intermediate )
  
  ## Histology
  # the his-function is a continuous version of histology (the unobserved truth)
  l <- rnorm(n, 0.5, 1) # unobserved risk factor
  his <- sqrt(x7*x6 + x1) + 2*l # ~15-20% of his is explained by l
  
  # let about share.lc% of all individuals get LC:
  cutoff <- quantile(his, probs = 1 - share.lc)
  lc.idx <- his >= cutoff
  g.lc   <- his[lc.idx]
  his[!lc.idx] <- 0 # no LC = zero histology
  
  ## create x8, which is what we are going to report (i.e. 4-level histology)
  # among the ones that got LC, let 50% be ADN, 25% SQM, 15% SCLC, 10% NSCLC.
  # let the order be (least to most harmful): ADN, SQM, NSCLC, SCLC
  borders   <- quantile(g.lc, probs = c(0.5, 0.75, 0.85))
  lc.histol <- rep(NA, length(g.lc))
  lc.histol[g.lc < borders[1]]                      <- 'ADN'
  lc.histol[borders[1] <= g.lc & g.lc < borders[2]] <- 'SQM'
  lc.histol[borders[2] <= g.lc & g.lc < borders[3]] <- 'NSCLC'
  lc.histol[g.lc >= borders[3]]                     <- 'SCLC' 
  x8          <- rep(NA, n)
  x8[lc.idx]  <- lc.histol
  x8[!lc.idx] <- '0' # won't develop LC
  
  
  ### 2) Nuisance parameter \nu ----
  # f is a function for the core factors of dying of lung cancer
  f <- 1/3 * x6 + 1/6*x7 + 0.2 * sqrt(x6*x7) + 0.5 * log(x1) 
  f <- (exp(f)^0.1) + sqrt(pi*2)
  
  # h is a function for the less important factors
  h          <- rep(NA, n)
  h[x3 <= 3] <- pi
  h[x3 > 3]  <- 2
  h          <- h + ifelse(x2 == 1, 2, 2.5) 
  nu         <- (h + f)^1.1 
  
  ### 3) Treatment effect tau ----
  # \tau is the treatment effect if Y is continuous. Otherwise, it is the driving factor behind the treatment effect if the outcome 
  
  # standardize x1, x2 and x7: in the sense that we need to map it in (0,1) space
  x1.st <- pnorm(x1, mean(x1), sqrt(var(x1)))
  x6.st <- pnorm(x6, mean(x6), sqrt(var(x6))) 
  x7.st <- pnorm(x7, mean(x7), sqrt(var(x7))) 
  
  # the following logistic functions require x \in (0,1) to return heterogeneous values
  logistic.fun.1 <- function(x) 1 / (1 + exp(-20 * (x - 1/3))) 
  logistic.fun.2 <- function(x) 2 / (1 + exp(-12 * (x - 1/2))) # sharper spike in the x = 1 region
  
  # continuous histology:
  # addition of noise term: there is an unobserved factor (which we have contained in the categorical histology)
  tau0        <- (1 + logistic.fun.2(x1.st)/3) * (1 + logistic.fun.2(x6.st)/3) * (1 + logistic.fun.1(x7.st)/2) + l/4 
  tau         <- rep(0, n)
  tau[lc.idx] <- - tau0[lc.idx]
  
  
  ### 4) Error term \varepsilon ----
  err <- rnorm(n, 0, 1)
  if(reporting.bias){
    sig.y <- 1/2 * sqrt(var(nu)) # standard deviation of error
    
    # now, assume that there is a reporting error: let there be 4 quantiles of x6 and each of them shares one error component. But: they get more intense with increasing values of x6.
    quants <- quantile(x6)[-1]
    err.group <- rep(NA, n)
    group <- rep(NA, n)
    for(j in 1:length(quants)){
      err.g <- abs(rnorm(1, j-1, 0.5))
      if(j == 1){
        err.group[x6 <= quants[j]] <- err.g
        group[x6 <= quants[j]] <- as.character(j)
      } else{
        err.group[quants[j-1] < x6 & x6 <= quants[j]] <- err.g
        group[quants[j-1] < x6 & x6 <= quants[j]] <- as.character(j)
      }
    }
    eps <- sig.y * (0.7*err + 0.3*err.group)
  } else{
    eps <- err
  }
  
  
  ### 5) Continuous (potential) outcomes ----
  # for the unlikely case of negative outcomes, set their lower bound to zero for the sake of logical consistency.
  y0.con <- nu + eps       ; y0.con[y0.con < 0] <- 0
  y1.con <- tau + nu + eps ; y1.con[y1.con < 0] <- 0
  w      <- rbinom(n, 1, 0.5)      # treatment status; assume RCT
  y.con  <- ifelse(w == 1, y1.con, y0.con) # observed outcome
  
  
  ### 6) Binary (potential) outcomes ---- 
  # moments of y0 -> only consider individuals that actually get LC! Such w/o LC shall have a zero chance of dying of LC throughout!
  mu    <- mean(y0.con[lc.idx])
  sigma <- sqrt(var(y0.con[lc.idx]))
  
  # standardize:
  y0.sd  <- (y0.con[lc.idx] - mu) / sigma
  y1.sd  <- (y1.con[lc.idx] - mu) / sigma
  
  # \pi = Pr(Y=1) for both cases
  pi0 <- pi1 <- rep(0, n)
  pi0[lc.idx] <- pnorm(y0.sd)
  pi1[lc.idx] <- pnorm(y1.sd)
  
  # create binary outcomes
  y0.bin <- rbinom(n, 1, pi0)
  y1.bin <- rbinom(n, 1, pi1)
  y.bin  <- ifelse(w == 1, y1.bin, y0.bin) # observed outcome
  
  
  
  ### 7) True parameter values ----
  
  #adjust life-years for lung cancer mortality here: assume 20% reduction in life-years left from lung cancer mortality
  Lifeyears = Lifeyears_intermediate *ifelse(y.bin == 1, 0.8, 1) 
  
  ate.con <- mean(y1.con) - mean(y0.con)
  hte.con <- y1.con - y0.con
  ate.bin <- mean(pi1) - mean(pi0)
  hte.bin <- pi1 - pi0
  rateratio.bin <-  rateratio.test(c(sum(y1.bin),sum(y0.bin)),c(sum(Lifeyears[w==1]),sum(Lifeyears[w==0]))) #NB: close to mean(pi1) / mean(pi0)
  
  ### 8) Organize Output in a List ----
  y        <- y.con ; y0 <- y0.con ; y1 <- y1.con # to be consistent with naming
  lst.con <- list(outcomes = data.frame(y, y0, y1), hte = hte.con, ate = ate.con)
  y       <- y.bin ; y0 <- y0.bin ; y1 <- y1.bin 
  lst.bin <- list(outcomes = data.frame(y, y0, y1), hte = hte.bin, ate = ate.bin, pi0 = pi0, pi1 = pi1,rateratio=rateratio.bin$estimate[1])
  
  data <- list(X = data.frame(x1, x2, x3, x4, x5, x6, x7, x8), 
               w = w, continuous.histology = his,
               continuous.outcomes = lst.con,
               binary.outcomes     = lst.bin)
  return(data)
}
