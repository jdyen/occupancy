#' @name occupancy
#'
#' @import jagsUI
#' 
#' @title fit occupancy-detection models
#'
#' @description occupancy is a function to fit occupancy-detection models in JAGS from
#'     within R. Models are specified with a formula interface and are supported by 
#'     methods to summarise, visualise, validate, and predict from fitted models.
#'
#' @param formula_occ model formula for occupancy component of the model. A two-sided
#'     formula with the response variable on the left of the \code{~} operator and the
#'     predictor variables on the right. Multiple predictor variables are separated by
#'     \code{+} operators. Random effects are included with vertical bars, using the
#'     notation \code{(1 | group)} to specify a random intercept for the variable
#'     \code{group}. More complex random structures (e.g., random slopes) are not 
#'     supported.
#'     
#' @param formula_detect model formula for detection component of the model. A one-sided
#'     formula with predictor variables on the right, formatted as for \code{formula_occ}.
#'     
#' @param site_id the name of the column in \code{data} in which site
#'     identifiers are recorded.
#'     
#' @param survey_id the name of the column in \code{data} in which survey
#'     identifiers are recorded.
#'     
#' @param data a \code{data.frame} in long format containing data on all variables. Required
#'     variables include the response (detection-nondetection data), site and survey identifiers
#'     (see \code{site_id} and \code{survey_id}, above), and predictor variables. Column names
#'     must match the names used in \code{formula_occ} and \code{formula_detect}.
#'     
#' @param jags_settings optional list of MCMC settings. Any or all items can be altered
#'     if needed. Options are:
#'     \describe{
#'       \item{\code{n_iter}}{the total number of MCMC iterations (including burn-in)}
#'       \item{\code{n_burnin}}{the number of MCMC iterations to discard as a burn-in}
#'       \item{\code{n_chains}}{the number of MCMC chains}
#'       \item{\code{n_thin}}{thinning rate of MCMC samples}
#'       \item{\code{parallel}}{logical, should chains be run in parallel?}
#'       \item{\code{modules}}{JAGS modules to load}
#'       \item{\code{params}}{character vector of parameters to store}
#'       \item{\code{seed}}{seed used to initialise MCMC chains}
#'     }
#'
#' @details This function fits an occupancy-detection model in JAGS from two formulas:
#'     \code{formula_occ} and \code{formula_detect}. Occupancy-detection models separate
#'     the two processes of being present at a site and being detect given presence at
#'     a site. This requires data from repeated visits (surveys) to sites.
#'     
#'     The occupancy component of the model (presence at a site) is defined at
#'     the site level. The detection component of the model (detections given presence)
#'     is defined at the survey level. The model assumes that associations between
#'     occupancy/detection and predictor variables are linear on a logit scale.
#'
#' @return \code{occupancy_model} - a \code{list} object that can be analysed using
#'   functions described in \link[occupancy:methods]{methods}.
#'   
#' @examples
#' \dontrun{
#'
#' # fit a model to simulated data
#' mod <- occupancy(response ~ occ_predictor1 + occ_predictor2 + 
#'                     (1 | occ_random1) + (1 | occ_random2),
#'                   ~ detect_predictor1 + detect_predictor2 + 
#'                    (1 | detect_random1),
#'                site_id = "site",
#'                survey_id = "survey",
#'                data = occupancy_data,
#'                jags_settings = list(n_iter = 1000, n_burnin = 500, n_thin = 2))
#'                
#' # plot the model coefficients
#' par(mfrow = c(2, 1))
#' plot(mod)
#' 
#' # extract the model coefficients
#' coef(mod)
#' 
#' # check model fit
#' calculate_metrics(mod)
#'      
#' }
#' 
#' @export
occupancy <- function(formula_occ, formula_detect, site_id, survey_id, data, jags_settings = list()) {
  
  # unpack jags_settings
  jags_set <- list(n_iter = 2000,
                   n_burnin = 1000,
                   n_chains = 4,
                   n_thin = 1,
                   parallel = FALSE,
                   modules = c("glm", "dic"),
                   params = c("beta_psi", "beta_p", "p_obs"),
                   seed = floor(runif(1, 1, 1e4)))
  jags_set[names(jags_settings)] <- jags_settings
  
  # unpack formulas and check variables exist
  data_tmp <- fetch_data(formula_occ, formula_detect, site_id, survey_id, data)
  var_names_occ <- data_tmp$var_names_occ
  var_names_detect <- data_tmp$var_names_detect
  data_tmp$var_names_occ <- NULL
  data_tmp$var_names_detect <- NULL
  
  # calculate required indices (nx_occ, nx_detect, nsite, nsurvey)
  nx_occ <- ncol(data_tmp$X_occ)
  nx_detect <- ncol(data_tmp$X_detect)
  nsite <- nrow(data_tmp$y)
  nsurvey <- ncol(data_tmp$y)
  nz_occ <- nz_detect <- ncluster_occ <- ncluster_detect <- NULL
  if (!is.null(data_tmp$Z_occ)) {
    nz_occ <- ncol(data_tmp$Z_occ)
    ncluster_occ <- apply(data_tmp$Z_occ, 2, max)
  }
  if (!is.null(data_tmp$Z_detect)) {
    nz_detect <- ncol(data_tmp$Z_detect)
    ncluster_detect <- apply(data_tmp$Z_detect, 2, max)
  }
  
  # check data dimensions are all OK
  with(data_tmp, check_dims(y, X_occ, X_detect, Z_occ, Z_detect))
  
  # which model?
  if (is.null(data_tmp$Z_occ) | is.null(data_tmp$Z_detect)) {
    if (is.null(data_tmp$Z_occ) & !is.null(data_tmp$Z_detect))
      jags_string <- detect_random()
    if (!is.null(data_tmp$Z_occ) & is.null(data_tmp$Z_detect))
      jags_string <- occ_random()
    if (is.null(data_tmp$Z_occ) & is.null(data_tmp$Z_detect))
      jags_string <- no_random()
  } else {
    jags_string <- all_random()
  }
  jags_set$model_file <- jags_string 

  # collate JAGS data set
  jags_data <- c(data_tmp,
                 list(nx_occ = nx_occ,
                      nx_detect = nx_detect,
                      nz_occ = nz_occ,
                      ncluster_occ = ncluster_occ,
                      nz_detect = nz_detect,
                      ncluster_detect = ncluster_detect,
                      nsite = nsite,
                      nsurvey = nsurvey))
  
  # remove missing data types
  jags_data <- jags_data[sapply(jags_data, function(x) !is.null(x))]
  
  # add some extra params if random effects are included
  if (!is.null(data_tmp$Z_occ))
    jags_set$params <- c(jags_set$params, "sigma_occ")
  if (!is.null(data_tmp$Z_detect))
    jags_set$params <- c(jags_set$params, "sigma_detect")
  
  # set some initial values
  inits <- function() {
    list(beta_psi = rep(0, jags_data$nx_occ),
         beta_p = rep(0, jags_data$nx_detect),
         z = apply(jags_data$y, 1, max))
  }

  # compile and run JAGS model
  out <- jags(data = jags_data,
              inits = inits,
              parameters.to.save = jags_set$params,
              model.file = textConnection(jags_set$model_file),
              n.chains = jags_set$n_chains,
              n.iter = jags_set$n_iter,
              n.burnin = jags_set$n_burnin,
              n.thin = jags_set$n_thin,
              parallel = jags_set$parallel,
              seed = jags_set$seed)
  
  
  # add the data to the output to allow model fit estimates,
  #    plotting and predictions
  out <- list(jags_model = out,
              settings = jags_set,
              data = c(jags_data,
                       var_names_occ = list(var_names_occ),
                       var_names_detect = list(var_names_detect)))
  
  # set a class so we can have default functions on top of this
  class(out) <- c("occupancy_model", class(out))  

  # return outputs
  out

}
