#' ----------------------------------------------------------------------
#' Summarizes the simulations from script "sim-grf.R. 
#' Makes plots and a csv tables.
#' 
#' Author: mwelz
#' Last changed: Feb. 25, 2021
#' ----------------------------------------------------------------------
rm(list = ls()) ; cat("\014")

# load the helper functions
source(paste0(getwd(), "/funs/estimation-funs.R"))

# load results of the simulation (from script "sim-grf.R")
load(paste0(getwd(), "/simulations/baseline/sim-grf/sim-results-raw.Rdata"))


### 1. make csv table for sample-level statistics ----

# overall statistics
final.table <- colMeans(grf.mat)

final.table <- rbind(c("Predicted absolute benefit",
                   "|Error predicted absolute benefit|",
                   "Coverage at 95% level",
                   "C index"),
                 round(final.table, 4))

write.table(final.table, 
            file = paste0(getwd(), "/simulations/baseline/sim-grf/sim-results.csv"),
            col.names = FALSE, row.names = FALSE, sep = ",")


### 2. make plot for group-level statistics ----
# make plots
library(ggplot2)
library(ggpmisc)

quantile.groups <- names(grf.arr.ls)

for(group in quantile.groups){
  
  ### 2.1 risk model, absolute benefit ----
  grf.ave        <- t(apply(grf.arr.ls[[group]], c(1,2), mean))
  risk.quantile <- factor(rownames(grf.ave))
  lv            <- levels(risk.quantile)
  lv            <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  grf.df        <- data.frame(grf.ave, risk.quantile)
  grf.table.abs <- data.frame(risk.quantile = rownames(grf.df), coverage = grf.df$abs.cover)

  gg <- ggplot(mapping = aes(x = abs.pb,
                             y = abs.ob, color = risk.quantile), data = grf.df) +
    annotate(geom = "table", x = -0.4, y = 0.5, label = list(grf.table.abs), 
             vjust = 1, hjust = 0) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = c(-0.4, 0.5), ylim = c(0, 0.5)) +
    labs(x = "Predicted absolute benefit",
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(paste0("GRF with ",
                   length(quantile.groups.ls[[group]]) + 1,
                   " risk groups (absolute)")) +
    theme(legend.position = "bottom")
  
  
  # open pdf
  pdf(file = paste0(getwd(), "/simulations/baseline/sim-grf/plots/",
                    "sim-summary_",
                    length(quantile.groups.ls[[group]]) + 1,
                    "-groups_grf_absolute.pdf"))
  print(gg)
  dev.off() # close pdf
  
} # FOR


# # subgroup analyses
# grf.subgroup <- round(apply(grf.arr, c(1,2), mean), 5)
# grf.subgroup <- cbind(c("Causal forest: Predicted absolute benefit",
#                        "Causal forest: Observed absolute benefit"),
#                       grf.subgroup)
# 
# grf.subgroup <- rbind(colnames(grf.subgroup), grf.subgroup)
# grf.subgroup[1,1] <- "Quantile baseline risk"
# 
# # append to final table
# overall <- cbind(overall, matrix("", 2, 1))
# final.table <- rbind(overall, rep("", 5),
#                      grf.subgroup)
