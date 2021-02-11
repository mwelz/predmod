rm(list = ls()) ; cat("\014")

# load results of the simulation (from script "sim-grf.R")
load(paste0(getwd(), "/simulations/baseline/sim-grf/sim-results-raw.Rdata"))

# overall statistics
overall <- colMeans(grf.mat)

overall <- rbind(c("Predicted absolute benefit",
                   "|Error predicted absolute benefit|",
                   "Coverage at 95% level",
                   "C index"),
                 round(overall, 4))

# subgroup analyses
grf.subgroup <- round(apply(grf.arr, c(1,2), mean), 5)
grf.subgroup <- cbind(c("Causal forest: Predicted absolute benefit",
                       "Causal forest: Observed absolute benefit"),
                      grf.subgroup)

grf.subgroup <- rbind(colnames(grf.subgroup), grf.subgroup)
grf.subgroup[1,1] <- "Quantile baseline risk"

# append to final table
overall <- cbind(overall, matrix("", 2, 1))
final.table <- rbind(overall, rep("", 5),
                     grf.subgroup)

write.table(final.table, 
            file = paste0(getwd(), "/simulations/baseline/sim-grf/sim-results.csv"),
            col.names = FALSE, row.names = FALSE, sep = ",")
