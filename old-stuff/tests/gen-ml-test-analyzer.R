rm(list = ls()) ; cat("\014")
library(ggplot2)

source(paste0(getwd(), "/funs/generic-ml/generic-ml-auxiliary-funs.R"))
load(file = paste0(getwd(), "/tests/gen-ml-test.Rdata"))

# there shouldn't be treatment effect heterogeneity!

# analyze
genML$VEIN$best.learners$GATES # difference is insignificant, so no hetero
genML$VEIN$best.learners$BLP  # beta2 is insignificant, so no hetero
genML$VEIN$best.learners$CLAN$z1 # there seems to be hetero along z1

## plot
# since treatment effect is negative, G1 is the most affected group (most negative), and G.K the least affected group (least negative)

# GATES
genericML.plot(genML, type = "GATES") # no hetero
genericML.plot(genML, type = "BLP")   # no hetero

genericML.plot(genML, type = "CLAN", CLAN.variable = "z1")   # no hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z2")   # hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z3")   # slight hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z4")   # no hetero
genericML.plot(genML, type = "CLAN", CLAN.variable = "z5")   # no hetero
