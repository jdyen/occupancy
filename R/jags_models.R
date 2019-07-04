all_random <- function() {
  "
model {
  
  # set priors for predictor matrices to account for missing data
  for (i in 1:nsite) {
    for (j in 1:nx_occ) {
      X_occ[i, j] ~ dnorm(0, 1)
    }
    for (j in 1:nx_detect) {
      for (k in 1:nsurvey) {
        X_detect[i, j, k] ~ dnorm(0, 1)
      }
    }
  }
  
  # set priors for regression coefficients on occupancy (psi) and detection (p)
  for (i in 1:nx_occ) {
    beta_psi[i] ~ dnorm(0, 0.1)
  }
  for (i in 1:nx_detect) {
    beta_p[i] ~ dnorm(0, 0.1)
  }
  
  # set priors for random effects in occupancy (psi) and detection (p)
  for (i in 1:nz_occ) {
    sigma_occ[i] ~ dunif(0, 1)
    prec_occ[i] <- pow(sigma_occ[i], -2)
    for (j in 1:ncluster_occ[i]) {
      gamma_occ[j, i] ~ dnorm(0, prec_occ[i])
    }
  }
  for (i in 1:nz_detect) {
    sigma_detect[i] ~ dunif(0, 1)
    prec_detect[i] <- pow(sigma_detect[i], -2)
    for (j in 1:ncluster_detect[i]) {
      gamma_detect[j, i] ~ dnorm(0, prec_detect[i])
    }
  }
  
  # occupancy (psi) has one value per site and depends on a matrix of predictors X_occ (defined at site level)
  for (i in 1:nsite)	{
    lpsi_tmp[i] <- inprod(beta_psi, X_occ[i, ])
    for (j in 1:nz_occ) {
      lpsi_tmp2[i, j] <- lpsi_tmp[i] + gamma_occ[Z_occ[i, j], j]
    }
    lpsi[i] <- sum(lpsi_tmp2[i, 1:nz_occ])
    logit(psi[i]) <- lpsi[i]
    z[i] ~ dbern(psi[i])
  }
  
  # detection has one value per survey and depends on a matrix of predictors X_detect (defined at survey level)
  for (i in 1:nsite) {
    for (j in 1:nsurvey) {
      
      # probability of detection depends on predictors
      lp_tmp[i, j] <- inprod(beta_p, X_detect[i, , j])
      for (k in 1:nz_detect) {
        lp_tmp2[i, j, k] <- lp_tmp[i, j] + gamma_detect[Z_detect[i, k, j], k]
      }
      lp[i, j] <- sum(lp_tmp2[i, j, 1:nz_detect])
      logit(p[i, j]) <- lp[i, j]
      
      # actual detection depends on occupancy and probability of detection
      p_obs[i, j] <- z[i] * p[i, j]
      y[i, j] ~ dbern(p_obs[i, j])
      
    }
  }
  
}
"
}

detect_random <- function() {
  "
model {

  # set priors for predictor matrices to account for missing data
  for (i in 1:nsite) {
  for (j in 1:nx_occ) {
  X_occ[i, j] ~ dnorm(0, 1e-4)
  }
  for (j in 1:nx_detect) {
  for (k in 1:nsurvey) {
  X_detect[i, j, k] ~ dnorm(0, 1e-4)
  }
  }
  }
  
  # set priors for regression coefficients on occupancy (psi) and detection (p)
  for (i in 1:nx_occ) {
  beta_psi[i] ~ dnorm(0, 1e-4)
  }
  for (i in 1:nx_detect) {
  beta_p[i] ~ dnorm(0, 1e-4)
  }
  
  # set priors for random effects in occupancy (psi) and detection (p)
  for (i in 1:nz_detect) {
  sigma_detect[i] ~ dunif(1e-4, 1000)
  prec_detect[i] <- pow(sigma_detect[i], -2)
  for (j in 1:ncluster_detect[i]) {
  gamma_detect[j, i] ~ dnorm(0, prec_detect[i])
  }
  }
  
  # occupancy (psi) has one value per site and depends on a matrix of predictors X_occ (defined at site level)
  for (i in 1:nsite)	{
  lpsi[i] <- inprod(beta_psi, X_occ[i, ])
  logit(psi[i]) <- lpsi[i]
  z[i] ~ dbern(psi[i])
  }
  
  # detection has one value per survey and depends on a matrix of predictors X_detect (defined at survey level)
  for (i in 1:nsite) {
  for (j in 1:nsurvey) {
  
  # probability of detection depends on predictors
  lp_tmp[i, j] <- inprod(beta_p, X_detect[i, , j])
  for (k in 1:nz_detect) {
  lp_tmp2[i, j, k] <- lp_tmp[i, j] + gamma_detect[Z_detect[i, k, j], k]
  }
  lp[i, j] <- sum(lp_tmp2[i, j, 1:nz_detect])
  logit(p[i, j]) <- lp[i, j]
  
  # actual detection depends on occupancy and probability of detection
  p_obs[i, j] <- z[i] * p[i, j]
  y[i, j] ~ dbern(p_obs[i, j])
  
  }
  }
  
}
"
}

occ_random <- function() {
  "
model {

  # set priors for predictor matrices to account for missing data
  for (i in 1:nsite) {
  for (j in 1:nx_occ) {
  X_occ[i, j] ~ dnorm(0, 1e-4)
  }
  for (j in 1:nx_detect) {
  for (k in 1:nsurvey) {
  X_detect[i, j, k] ~ dnorm(0, 1e-4)
  }
  }
  }
  
  # set priors for regression coefficients on occupancy (psi) and detection (p)
  for (i in 1:nx_occ) {
  beta_psi[i] ~ dnorm(0, 1e-4)
  }
  for (i in 1:nx_detect) {
  beta_p[i] ~ dnorm(0, 1e-4)
  }
  
  # set priors for random effects in occupancy (psi) and detection (p)
  for (i in 1:nz_occ) {
  sigma_occ[i] ~ dunif(1e-4, 1000)
  prec_occ[i] <- pow(sigma_occ[i], -2)
  for (j in 1:ncluster_occ[i]) {
  gamma_occ[j, i] ~ dnorm(0, prec_occ[i])
  }
  }
  
  # occupancy (psi) has one value per site and depends on a matrix of predictors X_occ (defined at site level)
  for (i in 1:nsite)	{
  lpsi_tmp[i] <- inprod(beta_psi, X_occ[i, ])
  for (j in 1:nz_occ) {
  lpsi_tmp2[i, j] <- lpsi_tmp[i] + gamma_occ[Z_occ[i, j], j]
  }
  lpsi[i] <- sum(lpsi_tmp2[i, 1:nz_occ])
  logit(psi[i]) <- lpsi[i]
  z[i] ~ dbern(psi[i])
  }
  
  # detection has one value per survey and depends on a matrix of predictors X_detect (defined at survey level)
  for (i in 1:nsite) {
  for (j in 1:nsurvey) {
  
  # probability of detection depends on predictors
  lp[i, j] <- inprod(beta_p, X_detect[i, , j])
  logit(p[i, j]) <- lp[i, j]
  
  # actual detection depends on occupancy and probability of detection
  p_obs[i, j] <- z[i] * p[i, j]
  y[i, j] ~ dbern(p_obs[i, j])
  
  }
  }
  
}
"
}

no_random <- function() {
  "
model {
  
  # set priors for predictor matrices to account for missing data
  for (i in 1:nsite) {
  for (j in 1:nx_occ) {
  X_occ[i, j] ~ dnorm(0, 1)
  }
  for (j in 1:nx_detect) {
  for (k in 1:nsurvey) {
  X_detect[i, j, k] ~ dnorm(0, 1)
  }
  }
  }
  
  # set priors for regression coefficients on occupancy (psi) and detection (p)
  for (i in 1:nx_occ) {
  beta_psi[i] ~ dnorm(0, 1)
  }
  for (i in 1:nx_detect) {
  beta_p[i] ~ dnorm(0, 1)
  }
  
  # occupancy (psi) has one value per site and depends on a matrix of predictors X_occ (defined at site level)
  for (i in 1:nsite)	{
  lpsi[i] <- inprod(beta_psi, X_occ[i, ])
  logit(psi[i]) <- lpsi[i]
  z[i] ~ dbern(psi[i])
  }
  
  # detection has one value per survey and depends on a matrix of predictors X_detect (defined at survey level)
  for (i in 1:nsite) {
  for (j in 1:nsurvey) {
  
  # probability of detection depends on predictors
  lp[i, j] <- inprod(beta_p, X_detect[i, , j])
  logit(p[i, j]) <- lp[i, j]
  
  # actual detection depends on occupancy and probability of detection
  p_obs[i, j] <- z[i] * p[i, j]
  y[i, j] ~ dbern(p_obs[i, j])
  
  }
  }
  
}
"
}