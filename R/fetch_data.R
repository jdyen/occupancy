# internal function to check formulas and extract data ready for JAGS
fetch_data <- function(formula_occ, formula_detect, site_id, survey_id, data) {
  
  # set up to ignore NAs if needed
  current_na_action <- options("na.action")$na.action
  options(na.action = "na.pass")
  
  # extract variables and data for occupancy model
  response <- all.vars(formula_occ)[1]
  terms <- terms(formula_occ)
  random_idx <- (grep("\\|", attributes(terms)$term.labels))
  var_names <- all.vars(delete.response(terms))
  occ_var_list <- colnames(attributes(terms)$factors)
  if (length(random_idx)) {
    occ_var_list_fixed <- occ_var_list[-grep("\\|", occ_var_list)]
  } else {
    occ_var_list_fixed <- occ_var_list
  }
  
  # use correct variable names when random is missing
  if (length(random_idx)) {
    # check there are no interactions in the random terms
    if (length(grep('\\*', occ_var_list[random_idx]))) {
      stop("cannot include interactions in random effects; use separate (1 | random) 
           terms for each random variable", call. = FALSE)
    }
    fixed_occ <- var_names[-random_idx]
    random_occ <- var_names[random_idx]
    } else {
      fixed_occ <- var_names
      random_occ <- NULL
    }
  
  # create x, y, z objects to pass to default method
  y <- get(response, envir = as.environment(data), inherits = TRUE)
  if (length(fixed_occ)) {
    x_tmp <- mget(fixed_occ, envir = as.environment(data), inherits = TRUE)
  }
  if (length(random_occ)) {
    z_tmp <- mget(random_occ, envir = as.environment(data), inherits = TRUE)
    z_tmp <- lapply(z_tmp, function(x) as.integer(as.factor(x)))
  }
  
  # create model matrix
  if (length(fixed_occ)) {
    X_occ <- model.matrix(as.formula(paste0("~", paste(occ_var_list_fixed, collapse = " + "))), data = x_tmp)
  } else {
    X_occ <- matrix(1, nrow = length(y), ncol = 1)
    colnames(X_occ) <- "intercept"
  }
  if (length(random_occ)) {
    Z_occ <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(random_occ, collapse = " + "))),
                          data = z_tmp)
  } else {
    Z_occ <- NULL
  }
  
  # repeat for detection model
  terms <- terms(formula_detect)
  random_idy <- (grep("\\|", attributes(terms)$term.labels))
  var_detect <- all.vars(delete.response(terms))
  detect_var_list <- colnames(attributes(terms)$factors)
  if (length(random_idy)) {
    detect_var_list_fixed <- detect_var_list[-grep("\\|", detect_var_list)]
  } else {
    detect_var_list_fixed <- detect_var_list
  }
  
  # use correct variable names when random is missing
  if (length(random_idy)) {
    # check there are no interactions in the random terms
    if (length(grep('\\*', detect_var_list[random_idy]))) {
      stop("cannot include interactions in random effects; use separate (1 | random) 
           terms for each random variable", call. = FALSE)
    }
    fixed_detect <- var_detect[-random_idy]
    random_detect <- var_detect[random_idy]
    } else {
      fixed_detect <- var_detect
      random_detect <- NULL
    }
  
  # create x and z objects to pass to default method
  if (length(fixed_detect)) {
    x_tmp <- mget(fixed_detect, envir = as.environment(data), inherits = TRUE)
  }
  if (length(random_detect)) {
    z_tmp <- mget(random_detect, envir = as.environment(data), inherits = TRUE)
    z_tmp <- lapply(z_tmp, function(x) as.integer(as.factor(x)))
  }
  
  # create model matrix
  if (length(fixed_detect)) {
    X_detect <- model.matrix(as.formula(paste0("~", paste(detect_var_list_fixed, collapse = " + "))), data = x_tmp)
  } else {
    X_detect <- matrix(1, nrow = length(y), ncol = 1)
    colnames(X_detect) <- "intercept"
  }
  if (length(random_detect)) {
    Z_detect <- model.matrix(as.formula(paste0(" ~ -1 + ", paste(random_detect, collapse = " + "))),
                             data = z_tmp)
  } else {
    Z_detect <- NULL
  }
  
  # pull out survey IDs so we can set up appropriate dims for each variable
  site_id <- get(site_id, envir = as.environment(data), inherits = TRUE)
  survey_id <- get(survey_id, envir = as.environment(data), inherits = TRUE)
  
  # work out required dims first
  nsite <- length(unique(site_id))
  survey_names <- unique(survey_id)
  nsurvey <- length(survey_names)
  
  # fix up y
  y_new <- tapply(y, list(site_id, survey_id), sum, default = NA)
  
  # now do occupancy predictors
  X_occ_new <- matrix(NA, nsite, ncol(X_occ))
  for (i in seq_len(ncol(X_occ)))
    X_occ_new[, i] <- tapply(X_occ[, i], site_id, mean_fn, na.rm = TRUE)
  if (!is.null(Z_occ)) {
    Z_occ_new <- matrix(NA, nsite, ncol(Z_occ))
    for (i in seq_len(ncol(Z_occ)))
      Z_occ_new[, i] <- tapply(Z_occ[, i], site_id, mean_fn, na.rm = TRUE)
  } else {
    Z_occ_new <- NULL
  }
  
  # and now detection predictors
  X_detect_new <- array(NA, dim = c(nsite, ncol(X_detect), nsurvey))
  for (i in seq_len(ncol(X_detect)))
    X_detect_new[, i, ] <- tapply(X_detect[, i], list(site_id, survey_id), mean_fn, na.rm = TRUE,
                                  simplify = TRUE)
  if (!is.null(Z_detect)) {
    Z_detect_new <- array(NA, dim = c(nsite, ncol(Z_detect), nsurvey))
    for (i in seq_len(ncol(Z_detect)))
      Z_detect_new[, i, ] <- tapply(Z_detect[, i], list(site_id, survey_id), mean_fn, na.rm = TRUE,
                                    simplify = TRUE, default = 1)
  } else {
    Z_detect_new <- NULL  
  } 
  
  # are there some occupancy predictors (defined at site level) that have more than one value per site?
  nval_occ <- matrix(NA, nsite, ncol(X_occ))
  for (i in seq_len(ncol(X_occ)))
    nval_occ[, i] <- tapply(X_occ[, i], site_id, function(x) length(unique(x[!is.na(x)], na.rm = TRUE)))
  nval_by_col <- apply(nval_occ, 2, function(x) any(x > 1))
  if (any(nval_by_col)) {
    col_ids <- colnames(X_occ)[nval_by_col]
    if (length(col_ids) > 1)
      col_ids <- paste(col_ids, collapse = ", ")
    warning("Some sites have multiple values for ", col_ids, call. = FALSE)
  }
  
  # are any random variables missing intermediate levels?
  if (!is.null(Z_occ_new))
    Z_occ_new <- apply(Z_occ_new, 2, function(x) as.integer(as.factor(x)))
  if (!is.null(Z_detect_new))
    Z_detect_new <- apply(Z_detect_new, c(2, 3), function(x) as.integer(as.factor(x)))
  
  # reset NA action
  options(na.action = current_na_action)
  
  # return outputs
  list(y = y_new,
       X_occ = X_occ_new, X_detect = X_detect_new,
       Z_occ = Z_occ_new, Z_detect = Z_detect_new,
       var_names_occ = occ_var_list_fixed, var_names_detect = detect_var_list_fixed)
  
}

# internal function to check dimensions of all data
check_dims <- function(y, X_occ, X_detect, Z_occ, Z_detect) {
  
  # check one row of occupancy predictors per site
  if (!is.null(X_occ)) {
    if (nrow(X_occ) != nrow(y)) {
      stop("Occupancy predictors should have one row per site but there are ",
           nrow(y), "sites and ", nrow(X_occ), " rows of occupancy predictors",
           call. = FALSE)
    }
  }
  
  # check one row of detection predictors per site (multiple slices for surveys)
  if (!is.null(X_detect)) {
    if (nrow(X_detect) != nrow(y)) {
      stop("Detection predictors should have one row per site but there are ",
           nrow(y), "sites and ", nrow(X_detect), " rows of detection predictors",
           call. = FALSE)
    }
  } 
  
  # check one slice of detection predictors per survey
  if (!is.null(X_detect)) {
    if (dim(X_detect)[3] != ncol(y)) {
      stop("There should be one set of detection predictors for each survey",
           call. = FALSE)
    } 
  }
  
  # check one row of occupancy random effects per site
  if (!is.null(Z_occ)) {
    if (nrow(Z_occ) != nrow(y)) {
      stop("Occupancy random effects should have one row per site but there are ",
           nrow(y), " sites and ", nrow(Z_occ), " rows of occupancy random effects",
           call. = FALSE)
    }
  } 
  
  # check one row of detection random effects per site (multiple slices for surveys)
  if (!is.null(Z_detect)) {
    if (nrow(Z_detect) != nrow(y)) {
      stop("Detection random effects should have one row per site but there are ",
           nrow(y), " sites and ", nrow(Z_detect), " rows of detection random effects",
           call. = FALSE)
    }
  } 
  
  # check one slice of detection random effects per survey
  if (!is.null(Z_detect)) {
    if (dim(Z_detect)[3] != ncol(y)) {
      stop("There should be one set of detection random effects for each survey",
           call. = FALSE)
    } 
  }
  
}
