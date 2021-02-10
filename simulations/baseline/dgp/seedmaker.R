rm(list = ls()) ; cat("\014")
R <- 100

set.seed(83267)
seeds <- round(1e7 * runif(R), 0)

save(seeds, 
     file = paste0(getwd(), "/simulations/baseline/dgp/seeds.Rdata"))