# This example code computes APS for a cutoff rule on three variables.
# The import requires ivaps.R to be placed in the same folder as this file.

# Edit this file as needed.
library(haven)
source("ivaps.R")

predict <- function(X) {
  return(as.integer((X[, 1] >= 0.202) & (X[, 2] >= 25000) & (X[,3] <= 0.03)))
}

df <- read_dta("data/safety_net_elig.dta")
df$aps <- estimate_aps(predict = predict, X = df[c("sum_pctg_ssi_mdcd_days", "ucc_per_bed", "profit_margin")], C = c(1,2,3), S = 1000, delta = 0.05, nprocesses = 2)
write_dta(df, "data/safety_net_elig.dta")
