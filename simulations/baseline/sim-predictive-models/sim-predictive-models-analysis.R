rm(list = ls()) ; cat("\014")

# load results of the simulation (from script "sim-predictive-models.R")
load(paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results-raw.Rdata"))

# overall statistics
overall <- rbind(colMeans(rm.mat), colMeans(em.mat))
rownames(overall) <- c("Risk model", "Effect model")

overall <- rbind(c("Predicted absolute benefit",
                   "|Error predicted absolute benefit|",
                   "Predicted relative benefit",
                   "|Error predicted relative benefit|",
                   "C index"),
                 round(overall, 5))
overall <- cbind(c("", "Risk model", "Effect model"), overall)

# prepare the final table
final.table <- cbind(overall, matrix("", 3, 6))
final.table <- rbind(final.table, matrix("", 2, ncol(final.table)))
final.table[5,1] <- "Frequency of selected variables in effect model"

# frequency of the selected variables by the effect model
em.varselection <- colMeans(em.vars)
final.table <- rbind(final.table, c("Variable in effect model", 
                                    names(em.varselection)))
final.table <- rbind(final.table, c("Frequency", em.varselection))


# subgroup analyses
rm.subgroup <- round(apply(rm.arr, c(1,2), mean), 5)
rm.subgroup <- cbind(c("Risk model: Predicted absolute benefit",
                       "Risk model: Observed absolute benefit",
                       "Risk model: Predicted relative benefit",
                       "Risk model: Observed relative benefit"), 
                     rm.subgroup)
rm.subgroup <- rbind(colnames(rm.subgroup), rm.subgroup)
rm.subgroup[1,1] <- "Quantile baseline risk"

em.subgroup <- round(apply(em.arr, c(1,2), mean), 5)
em.subgroup <- cbind(c("Effect model: Predicted absolute benefit",
                       "Effect model: Observed absolute benefit",
                       "Effect model: Predicted relative benefit",
                       "Effect model: Observed relative benefit"), 
                     em.subgroup)
em.subgroup <- rbind(colnames(em.subgroup), em.subgroup)
em.subgroup[1,1] <- "Quantile baseline risk"

# Combine subgroups
subgroups <- rbind(rm.subgroup, rep("", ncol(rm.subgroup)), em.subgroup)
subgroups <- cbind(subgroups, matrix("", nrow(subgroups), 7))

# append to final table
final.table <- rbind(final.table, 
                     rep("", ncol(final.table)), rep("", ncol(final.table)), 
                     subgroups)



write.table(final.table, 
            file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results.csv"),
            col.names = FALSE, row.names = FALSE, sep = ",")

