#' ----------------------------------------------------------------------
#' Summarizes the simulations from script "sim-predictive-models.R. 
#' Makes plots and a csv tables.
#' 
#' Author: mwelz
#' Last changed: Feb. 25, 2021
#' ----------------------------------------------------------------------
rm(list = ls()) ; cat("\014")

# load the helper functions 
source(paste0(getwd(), "/funs/estimation-funs.R"))

# load results of the simulation (from script "sim-predictive-models.R")
load(paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results-raw.Rdata"))

### 1. make csv table for sample-level statistics ----

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

write.table(final.table,
            file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results.csv"),
            col.names = FALSE, row.names = FALSE, sep = ",")



### 2. make plot for group-level statistics ----
# make plots
library(ggplot2)
library(ggpmisc)

quantile.groups <- names(rm.arr.ls)

for(group in quantile.groups){
  
  ### 2.1 risk model, absolute benefit ----
  rm.ave        <- t(apply(rm.arr.ls[[group]], c(1,2), mean))
  risk.quantile <- factor(rownames(rm.ave))
  lv            <- levels(risk.quantile)
  lv            <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  rm.df        <- data.frame(rm.ave, risk.quantile)
  rm.table.abs <- data.frame(risk.quantile = rownames(rm.df), coverage = rm.df$abs.cover)
  rm.table.rel <- data.frame(risk.quantile = rownames(rm.df), coverage = rm.df$rel.cover)
  
  gg <- ggplot(mapping = aes(x = abs.pb,
                             y = abs.ob, color = risk.quantile), data = rm.df) +
    annotate(geom = "table", x = -0.4, y = 0.5, label = list(rm.table.abs), 
             vjust = 1, hjust = 0) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = c(-0.4, 0.5), ylim = c(0, 0.5)) +
    labs(x = "Predicted absolute benefit",
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(paste0("Risk model with ",
                   length(quantile.groups.ls[[group]]) + 1,
                   " risk groups (absolute)")) +
    theme(legend.position = "bottom")
  
  
  # open pdf
  pdf(file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/plots/",
                    "sim-summary_",
                    length(quantile.groups.ls[[group]]) + 1,
                    "-groups_risk-model_absolute.pdf"))
  print(gg)
  dev.off() # close pdf
    
  ### 2.2 risk model, relative benefit ----
  
  gg <- ggplot(mapping = aes(x = rel.pb,
                             y = rel.ob, color = risk.quantile), data = rm.df) +
    annotate(geom = "table", x = 0.8, y = 0.8, label = list(rm.table.rel), 
             vjust = 1, hjust = 0) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = c(0.5, 1), ylim = c(0.5, 0.8)) +
    labs(x = "Predicted relative benefit",
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(paste0("Risk model with ",
                   length(quantile.groups.ls[[group]]) + 1,
                   " risk groups (relative)")) +
    theme(legend.position = "bottom")
  
  # open pdf
  pdf(file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/plots/",
                    "sim-summary_",
                    length(quantile.groups.ls[[group]]) + 1,
                    "-groups_risk-model_relative.pdf"))
  print(gg)
  dev.off() # close pdf
  
  
  ### 2.1 effect model, absolute benefit ----
  em.ave        <- t(apply(em.arr.ls[[group]], c(1,2), mean))
  risk.quantile <- factor(rownames(em.ave))
  lv            <- levels(risk.quantile)
  lv            <- lv[order.intervals(lv, quantile.nam = TRUE)]
  risk.quantile <- factor(risk.quantile, levels = lv)
  
  em.df        <- data.frame(em.ave, risk.quantile)
  em.table.abs <- data.frame(risk.quantile = rownames(em.df), coverage = em.df$abs.cover)
  em.table.rel <- data.frame(risk.quantile = rownames(em.df), coverage = em.df$rel.cover)
  
  gg <- ggplot(mapping = aes(x = abs.pb,
                             y = abs.ob, color = risk.quantile), data = em.df) +
    annotate(geom = "table", x = -0.4, y = 0.5, label = list(em.table.abs), 
             vjust = 1, hjust = 0) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = c(-0.4, 0.5), ylim = c(0, 0.5)) +
    labs(x = "Predicted absolute benefit",
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(paste0("Effect model with ",
                   length(quantile.groups.ls[[group]]) + 1,
                   " risk groups (absolute)")) +
    theme(legend.position = "bottom")
  
  # open pdf
  pdf(file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/plots/",
                    "sim-summary_",
                    length(quantile.groups.ls[[group]]) + 1,
                    "-groups_effect-model_absolute.pdf"))
  print(gg)
  dev.off() # close pdf
  
  ### 2.2 effect model, relative benefit ----

  gg <- ggplot(mapping = aes(x = rel.pb,
                             y = rel.ob, color = risk.quantile), data = em.df) +
    annotate(geom = "table", x = 0.8, y = 0.8, label = list(em.table.rel), 
             vjust = 1, hjust = 0) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = c(0.5, 1), ylim = c(0.5, 0.8)) +
    labs(x = "Predicted relative benefit",
         y = "Observed benefit") +
    theme_classic() +
    ggtitle(paste0("Effect model with ",
                   length(quantile.groups.ls[[group]]) + 1,
                   " risk groups (relative)")) +
    theme(legend.position = "bottom")
  
  # open pdf
  pdf(file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/plots/",
                    "sim-summary_",
                    length(quantile.groups.ls[[group]]) + 1,
                    "-groups_effect-model_relative.pdf"))
  print(gg)
  dev.off() # close pdf
    
} # END for


# # subgroup analyses
# rm.subgroup <- round(apply(rm.arr, c(1,2), mean), 5)
# rm.subgroup <- cbind(c("Risk model: Predicted absolute benefit",
#                        "Risk model: Observed absolute benefit",
#                        "Risk model: Predicted relative benefit",
#                        "Risk model: Observed relative benefit"), 
#                      rm.subgroup)
# rm.subgroup <- rbind(colnames(rm.subgroup), rm.subgroup)
# rm.subgroup[1,1] <- "Quantile baseline risk"
# 
# em.subgroup <- round(apply(em.arr, c(1,2), mean), 5)
# em.subgroup <- cbind(c("Effect model: Predicted absolute benefit",
#                        "Effect model: Observed absolute benefit",
#                        "Effect model: Predicted relative benefit",
#                        "Effect model: Observed relative benefit"), 
#                      em.subgroup)
# em.subgroup <- rbind(colnames(em.subgroup), em.subgroup)
# em.subgroup[1,1] <- "Quantile baseline risk"
# 
# # Combine subgroups
# subgroups <- rbind(rm.subgroup, rep("", ncol(rm.subgroup)), em.subgroup)
# subgroups <- cbind(subgroups, matrix("", nrow(subgroups), 7))
# 
# # append to final table
# final.table <- rbind(final.table, 
#                      rep("", ncol(final.table)), rep("", ncol(final.table)), 
#                      subgroups)
# 
# 
# 
# write.table(final.table, 
#             file = paste0(getwd(), "/simulations/baseline/sim-predictive-models/sim-results.csv"),
#             col.names = FALSE, row.names = FALSE, sep = ",")
# 
