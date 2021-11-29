library(parallel)
library(ivpack)
library(dplyr)
library(car)
library(lmtest)
library(sandwich)

estimate_aps <- function (predict, X, C, S = 100, delta = 0.1, seed = 0, nprocesses = 1) {
  # Parameters
  # -----------
  # predict: function
  #     Function taking a 2D design matrix and returning a vector of predictions
  # X: data frame
  #     2D design matrix
  # C: vector
  #     Integer column indices for continuous variables
  # S: integer, default: 100
  #     Number of draws for each APS estimation
  # delta: double, default: 0.1
  #     Radius of sampling ball
  # seed: integer, default: 0
  #     Seed for random number generator
  # nprocesses: integer, default: 1
  #     Number of processes used to parallelize APS estimation
  # Returns
  # -----------
  # vector
  #     Vector of estimated APS for each observation in sample
  
  set.seed(seed)
  X = as.matrix(X)
  X_c = X[, C]
  if (length(C) == 1) {
    dim(X_c) = c(length(X_c), 1)
  }
  c_std = apply(X_c, 2, sd)
  
  samples = mclapply(seq(1, S), estimate_aps_helper, delta, X_c, c_std, X, C, predict, mc.cores = nprocesses)
  return(Reduce("+", samples)/length(samples))
}

estimate_aps_helper <- function (i, delta, X_c, c_std, X, C, predict) {
  # Resample continuous features
  dev = runif(length(X_c), -delta, delta)
  dim(dev) = dim(X_c)
  X_c_s = X_c + t(c_std * t(dev))
  X[, C] = X_c_s
  return(predict(X))
}

estimate_treatment_effect <- function (data, aps, Y, Z, D, W = NULL, cov_type = "HC1", weights = NULL) {
  # Parameters
  # -----------
  # data: data frame
  #     Data frame with relevant variables
  # aps: string
  #     Variable name for estimated APS
  # Y: string
  #     Variable name for outcomes
  # Z: string
  #     Variable name for treatment recommendations
  # D: string
  #     Variable name for treatment assignments
  # W: vector, default: NULL
  #     Variable names for controls
  # cov_type: string, default: "HC1"
  #     Covariance type of IV2SLS.
  # weights: string, default: NULL
  #     Variable name for observation weights used in estimation
  # Returns
  # -----------
  # coeftest object
  #     Matrix containing regression results.
  
  obs = nrow(data)
  # Use only observations where aps is nondegenerate
  if (aps == "aps") {
    data = subset(data, (aps != 0) & (aps != 1))
  } else {
    data = subset(data, (eval(parse(text=aps)) != 0) & (eval(parse(text=aps)) != 1))
  }
  print(paste("We will fit on", nrow(data), "values out of", obs, "from the dataset for which the APS estimation is nondegenerate."))
  
  # Check for single non-degeneracy
  single = dim(unique(data[c(aps)]))[1] == 1
  
  if (!is.null(weights)) {
    weights = data[c(weights)][,1]
  }
  
  if (single) {
    results = ivreg(eval(paste(Y, "~", paste(c(D, aps, W), collapse = " + "), "-1", "|", paste(c(Z, aps, W), collapse = " + "))), data = data, weights = weights)
  } else {
    results = ivreg(eval(paste(Y, "~", paste(c(D, aps, W), collapse = " + "), "|", paste(c(Z, aps, W), collapse = " + "))), data = data, weights = weights)
  }
  
  return(coeftest(results, vcov = vcovHC, type = cov_type))
}

covariate_balance_test <- function (data, aps, X, Z, W = NULL, cov_type = "HC1") {
  # Parameters
  # -----------
  # data: data frame
  #     Data frame with relevant variables
  # aps: string
  #     Variable name for estimated APS
  # X: vector
  #     Variable names for balance covariates
  # Z: string
  #     Variable name for treatment recommendations
  # W: vector, default: NULL
  #     Variable names for controls
  # cov_type: string, default: "HC1"
  #     Covariance type of regression
  # Returns
  # -----------
  # list
  #     X = coeftest object for balance regressions
  #     joint = linearHypothesis.mlm object for joint hypothesis test 
  
  obs = nrow(data)
  # Use only observations where aps is nondegenerate
  if (aps == "aps") {
    data = subset(data, (aps != 0) & (aps != 1))
  } else {
    data = subset(data, (eval(parse(text=aps)) != 0) & (eval(parse(text=aps)) != 1))
  }
  print(paste("We will fit on", nrow(data), "values out of", obs, "from the dataset for which the APS estimation is nondegenerate."))
  
  # Check for single non-degeneracy
  single = dim(unique(data[c(aps)]))[1] == 1
  
  if (single) {
    results = lm(eval(paste(paste0("cbind","(",paste(X, collapse = ", "),")"), "~", paste(c(Z, aps, W), collapse = " + "), "-1")), data = data)
  } else {
    results = lm(eval(paste(paste0("cbind","(",paste(X, collapse = ", "),")"), "~", paste(c(Z, aps, W), collapse = " + "))), data = data)
  }
  
  return(list("X" = coeftest(results, vcov = vcovHC, type = cov_type), "joint" = linearHypothesis(results, hypothesis.matrix = c(paste0(Z," = 0")))))
}
