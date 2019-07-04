#' @name validate
#'
#' @importFrom pROC auc roc
#' 
#' @title perform cross validation on a fitted occupancy-detection model
#'
#' @description validate is a function to cross validate fitted occupancy-detection
#'     models. This function will automatically prepare folds and run full cross validation.
#'     Note: cross validation requires fitting models repeatedly and can take some time.
#'
#' @param object occupancy model fitted with \link[occupancy:occupancy]{occupancy::occupancy}.
#'     
#' @param n_cv number of cross validation folds. Defaults to 10.
#'
#' @details Cross validation is a useful method to test model reliability on new data. Many
#'     estimates of model fit are susceptible to overfitting. Cross validation attempts
#'     to avoid this risk by repeatedly partitioning the data and generating predictions
#'     for parts of the data not used during model fitting. This is an approximation of
#'     true model predictive capacity and still entails some risk of overfitting because the
#'     data come from the same statistical population and are not fully independent from the
#'     data used to fit models.
#'
#' @return a \code{list} object with cross-validated estimates of McFadden's r2 and AUC.
#'   
#' @examples
#' \dontrun{
#'
#' NULL
#' 
#' }
NULL

#' @export
validate <- function(object, n_cv = 10) {
  
  n_obs <- object$data$nsite
  cv_size <- floor(n_obs / n_cv)
  cv_vec <- sample(seq_len(n_obs), size = n_obs, replace = FALSE)
  cv_list <- vector("list", length = n_cv)
  for (i in seq_len(n_cv)) {
    if (i < n_cv) {
      start <- (i - 1) * cv_size + 1
      end <- i * cv_size
      cv_list[[i]] <- cv_vec[start:end]
    } else {
      start <- (i - 1) * cv_size + 1
      cv_list[[i]] <- cv_vec[start:length(cv_vec)]
    }
  }
  
  # prepare a list of data for each fold
  data_list <- lapply(cv_list, prepare_cv_data, object$data)
  
  # prepare a list of newdata for predictions
  newdata_list <- lapply(cv_list, prepare_cv_newdata, object$data)
  
  # pull out settings
  settings <- object$settings
  
  # fit to all folds and return in a list  
  fitted <- lapply(seq_len(n_cv), cv_inner, settings, data_list, newdata_list)
  fitted <- do.call(rbind, fitted)
  
  # pull out fitted and observed values (need to calculate log likelihoods)  
  observed <- object$data$y[cv_vec, ]
  
  # keep fitted values within (0, 1)
  fitted[fitted > 0.999] <- 0.999
  fitted[fitted < 0.001] <- 0.001
  
  # pull out intercept so we can get a null loglikelihood  
  intercept_prob <- mean(fitted, na.rm = TRUE)
  
  # calculate loglikelihoods
  loglik_full <- sum(log(fitted) * observed + log(1 - fitted) * (1 - observed), na.rm = TRUE)
  loglik_null <- log(intercept_prob) * sum(observed, na.rm = TRUE) + log(1 - intercept_prob) * sum(1 - observed, na.rm = TRUE)
  
  # calculate r2 value
  r2 <- 1 - loglik_full / loglik_null
  
  # AUC calc
  auc <- pROC::auc(pROC::roc(c(observed), c(fitted)))
  
  # return outputs
  list(r2_cv = r2, auc_cv = auc)
  
}

# internal function to prepare CV folds
prepare_cv_data <- function(idx, data) {
  
  # pull out data without fold i
  data <- list(y = data$y[-idx, ],
               X_occ = data$X_occ[-idx, ],
               X_detect = data$X_detect[-idx, , ],
               Z_occ = data$Z_occ[-idx, ],
               Z_detect = data$Z_detect[-idx, , ],
               nx_occ = data$nx_occ,
               nx_detect = data$nx_detect,
               nz_occ = data$nz_occ,
               nz_detect = data$nz_detect,
               nsurvey = data$nsurvey)
  
  # double check factor levels in random effects
  if (data$nz_occ == 1)
    data$Z_occ <- matrix(data$Z_occ, ncol = 1)
  if (data$nz_detect == 1) {
    dims <- dim(data$Z_detect)
    data$Z_detect <- array(data$Z_detect, dim = c(dims[1], 1, dims[2]))
  }
  data$Z_occ <- apply(data$Z_occ, 2, function(x) as.integer(as.factor(x)))
  data$Z_detect <- apply(data$Z_detect, c(2, 3), function(x) as.integer(as.factor(x)))
  
  # add in remaining indices
  data$nsite <- nrow(data$y)
  data$ncluster_occ <- apply(data$Z_occ, 2, max)
  data$ncluster_detect <- apply(data$Z_detect, 2, max)
  
  # return outputs
  data
  
}

# internal function to prepare CV folds
prepare_cv_newdata <- function(idx, data) {
  
  # pull out data with only fold i
  out <- list(X_occ = data$X_occ[idx, ],
              X_detect = data$X_detect[idx, , ])
  
  # replace NAs with zeros (set missing to mean)
  out$X_occ <- ifelse(is.na(out$X_occ), 0, out$X_occ)
  out$X_detect <- ifelse(is.na(out$X_detect), 0, out$X_detect)
  
  # return data set
  out
  
}

# internal function to fit one CV fold
cv_inner <- function(i, settings, data_list, newdata_list) {
  
  # set some initial values
  inits <- function() {
    list(beta_psi = rep(0, data_list[[i]]$nx_occ),
         beta_p = rep(0, data_list[[i]]$nx_detect),
         z = apply(data_list[[i]]$y, 1, max))
  }
  
  # compile and run JAGS model
  mod <- jags(data = data_list[[i]],
              inits = inits,
              parameters.to.save = settings$params,
              model.file = settings$file_path,
              n.chains = settings$n_chains,
              n.iter = settings$n_iter,
              n.burnin = settings$n_burnin,
              n.thin = settings$n_thin,
              parallel = settings$parallel,
              seed = settings$seed)
  
  # generate predictions
  full_mod <- list(jags_model = mod,
                   data = c(data_list[[i]]))
  
  # set a class so we can have default functions on top of this
  class(full_mod) <- c("occupancy_model", class(full_mod))
  
  # make and return predictions
  predict(full_mod, newdata_list[[i]], type = "response")$predicted
  
}
