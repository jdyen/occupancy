#' @name methods
#'
#' @title methods for occupancy_models
#'
#' @importFrom graphics axis lines plot points
#' @importFrom stats as.formula coef delete.response model.matrix plogis predict quantile rbinom runif sd
#' @importFrom utils stack
#' @importFrom ggplot2 alpha
#'
#' @description This is a list of functions (mostly from base R) that are
#'   currently implemented for fitted occupancy models.
#'
#' @section Usage: \preformatted{
#'
#'  # predict if newdata are provided as a data.frame
#'  predict(object, newdata = NULL, type = c("link", "response"), ...)
#'  
#'  # predict if newdata are provided as a raster
#'  spatial_predict(object, newdata = NULL, type = "response", ...)
#'  
#'  # extract coefficients from a fitted model
#'  coef(object, ...)
#'  
#'  # summarise a fitted model
#'  summary(object, ...)
#'  
#'  # extract fitted values from a model
#'  fitted(object, type = c("link", "response"), ...)
#'  
#'  # calculate a pseudo-r2 value based on McFadden's r-squared
#'  r2_calc(object, ...)
#'  
#'  # calculate a suite of validation metrics (r2, AUC, DIC)
#'  calculate_metrics(object, ...)
#'  
#'  # plot the model coefficients
#'  plot(x, type, names_occ, names_detect, intercept, ...)
#'  
#'  # predictive plots of probabilities of occupancy
#'  plot_pr_occ(object, npred, var_name, label, ...)
#'  
#'  # predictive plots of probabilities of detection
#'  plot_pr_detect(object, npred, var_name, label, scale, ...)
#'  
#' }
#' 
#' @param object
#' @param newdata
#' @param type
#' @param \dots
#' @param x
#' @param names_occ
#' @param names_detect
#' @param intercept
#' @param npred
#' @param var_name
#' @param label
#' @param scale
#'
#' @details \code{predict} generates predictions of occupancy probabilities, detection
#'     probabilities, and likely detections (sampled as binary detection/nondetection)
#'     given these probabilities.
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
#' # spatial predictions
#' 
#' # plot_pr_occ example
#' 
#' # plot_pr_detect
#' 
#' }
NULL

#' @export
#' @rdname methods
#' 
predict.occupancy_model <- function(object, newdata = NULL, type = c("link", "response"), ...) {
  
  # check if link or response or neither (default to link)
  if (length(type) == 2)
    type <- "link"
  
  # pull out regression coefs
  betas <- lapply(object$jags_model$samples, function(x) x[, grep("beta", colnames(x))])
  betas <- do.call(rbind, betas)
  beta_p <- betas[, grep("beta_p\\[", colnames(betas))]
  beta_psi <- betas[, grep("beta_psi\\[", colnames(betas))]
  
  # predict fitted values if new data are not provided
  if (is.null(newdata))
    newdata <- list(X_occ = object$data$X_occ, X_detect = object$data$X_detect)
  
  # make predictions of occupancy
  out <- apply(beta_psi %*% t(newdata$X_occ), 2, mean)
  p_occ <- plogis(out)
  
  # is it actually present?
  out <- rbinom(length(out), prob = p_occ, size = 1)
  out <- matrix(out, nrow = nrow(newdata$X_occ))
  
  # now need to consider predictions of detectability
  out2 <- matrix(NA, nrow = nrow(out), ncol = ncol(object$data$y))
  for (i in seq_len(dim(newdata$X_detect)[3])) {
    out2[, i] <- apply(beta_p %*% t(newdata$X_detect[, , i]), 2, mean)
  }
  p_detect <- plogis(out2)
  
  # combine the two
  out <- sweep(out2, 1, out, "*")
  
  # convert to response scale if required
  if (type == "response")
    out <- plogis(out)
  
  # return outputs
  list(predicted = out,
       p_occ = p_occ,
       p_detect = p_detect)
  
}

#' @export
#' @rdname methods
#' 
spatial_predict <- function(object, newdata = NULL, type = "response", ...) {
  
  # what if newdata aren't provided?
  if (is.null(newdata))
    stop("spatial_predict requires raster inputs; use predict to calculate fitted values", call. = FALSE)
  
  # check the inputs are appropriate raster objects
  if (!class(newdata) %in% c("RasterStack", "RasterBrick"))
    stop("newdata must be a RasterStack or RasterBrick object")
  
  # pull out coefficients
  betas <- coef(object)$occupancy[, "Mean"]
  
  # are there enough layers in the raster data?
  if (nlayers(newdata) != (length(betas) - 1))
    stop("spatial_predict requires one layer per predictor variable")
  
  # add an intercept to the raster stack
  newdata <- stack(setValues(newdata[[1]], values = 1), newdata)
  
  # calculate effect of each variable and sum over all variables
  out <- setValues(newdata[[1]], values = 0)
  for (i in seq_along(betas))
    out <- out + betas[i] * newdata[[i]]
  
  # convert back to response (probability) scale if needed
  if (type == "response")
    out <- calc(out, plogis)
  
  # return outputs
  out
  
}

#' @export
#' @rdname methods
#' 
summary.occupancy_model <- function(object, ...) {
  
  summary(object$jags_model$samples, ...)
  
}

#' @export
fitted.occupancy_model <- function(object, type = c("link", "response"), ...) {
  
  if (type == "response") {
    mod_sum <- summary(object)
    fitted <- matrix(mod_sum$statistics[grep("p_obs", rownames(mod_sum$statistics)), "Mean"], ncol = 2)
  } else {
    fitted <- predict(object, type = "link")$predicted
  }
  
  fitted
  
}

#' @export
#' @rdname methods
#' 
coef.occupancy_model <- function(object, ...) {
  
  beta_psi <- lapply(object$jags_model$samples, function(y) y[, grep("beta_psi\\[", colnames(y))])
  beta_p <- lapply(object$jags_model$samples, function(y) y[, grep("beta_p\\[", colnames(y))])
  beta_psi <- do.call(rbind, beta_psi)
  beta_p <- do.call(rbind, beta_p)
  
  quantiles_psi <- t(apply(beta_psi, 2, quantile, p = c(0.025, 0.1, 0.5, 0.9, 0.975)))
  quantiles_psi <- cbind(apply(beta_psi, 2, mean), apply(beta_psi, 2, sd), quantiles_psi)
  quantiles_p <- t(apply(beta_p, 2, quantile, p = c(0.025, 0.1, 0.5, 0.9, 0.975)))
  quantiles_p <- cbind(apply(beta_p, 2, mean), apply(beta_p, 2, sd), quantiles_p)
  
  colnames(quantiles_p) <- colnames(quantiles_psi) <-
    c("Mean", "SD", "2.5%", "10%", "50%", "75%", "97.5%")
  
  rownames(quantiles_psi) <- c("Intercept", object$data$var_names_occ)
  rownames(quantiles_p) <- c("Intercept", object$data$var_names_detect)
  
  list(occupancy = quantiles_psi, detection = quantiles_p)
  
}

#' @export
#' @rdname methods
#' 
calculate_metrics <- function(object, ...) {
  
  # pull out fitted and observed values (need to calculate log likelihoods)  
  fitted <- fitted(object, type = "response")
  observed <- object$data$y
  
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
  
  # AIC/DIC calc?
  dic <- object$jags_model$DIC
  
  # return outputs
  list(r2 = r2, auc = auc, dic = dic)
  
}

#' @export
#' @rdname methods
#' 
r2_calc <- function(object, ...) {
  
  # pull out fitted and observed values (need to calculate log likelihoods)  
  fitted <- fitted(object, type = "response")
  observed <- object$data$y
  
  # keep fitted values within (0, 1)
  fitted[fitted > 0.999] <- 0.999
  fitted[fitted < 0.001] <- 0.001
  
  # pull out intercept so we can get a null loglikelihood  
  intercept_prob <- mean(fitted, na.rm = TRUE)
  
  # calculate loglikelihoods
  loglik_full <- sum(log(fitted) * observed + log(1 - fitted) * (1 - observed), na.rm = TRUE)
  loglik_null <- log(intercept_prob) * sum(observed, na.rm = TRUE) + log(1 - intercept_prob) * sum(1 - observed, na.rm = TRUE)
  
  # return r2 value
  1 - loglik_full / loglik_null
  
}

# internal function to support barplots
barplot_inner <- function(x, var_names, intercept, ...) {
  
  if (intercept)
    c("Intercept", var_names)
  
  if (!intercept)
    x <- x[, -1]
  
  nplot <- ncol(x)
  ylims <- range(x)
  xlims <- c(0.5, nplot + 0.5)
  plot(x[3, ] ~ seq_len(nplot),
       type = "n", bty = "l", xaxt = "n", yaxt = "n",
       xlim = xlims, ylim = ylims,
       xlab = "Variables", ylab = "Estimate", ...)
  for (i in seq_len(nplot)) {
    lines(c(i, i), c(x[1, i], x[5, i]), lwd = 2)
    lines(c(i, i), c(x[2, i], x[4, i]), lwd = 4)
  }
  points(x[3, ] ~ seq_len(nplot), pch = 16)
  lines(c(0, nplot + 1), c(0, 0), lty = 2, lwd = 2)
  axis(1, at = seq_len(nplot), labels = var_names, ...)
  axis(2, las = 1)
  
}

#' @export
#' @rdname methods
#' 
plot.occupancy_model <- function(x, type = "barplot", 
                                 names_occ = NULL, names_detect = NULL,
                                 intercept = FALSE,
                                 ...) {
  
  beta_psi <- lapply(x$jags_model$samples, function(y) y[, grep("beta_psi\\[", colnames(y))])
  beta_p <- lapply(x$jags_model$samples, function(y) y[, grep("beta_p\\[", colnames(y))])
  beta_psi <- do.call(rbind, beta_psi)
  beta_p <- do.call(rbind, beta_p)
  
  quantiles_psi <- apply(beta_psi, 2, quantile, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
  quantiles_p <- apply(beta_p, 2, quantile, p = c(0.025, 0.1, 0.5, 0.9, 0.975))
  
  if (is.null(names_occ))
    names_occ <- x$data$var_names_occ
  if (is.null(names_detect))     
    names_detect <- x$data$var_names_detect
  
  barplot_inner(quantiles_psi, names_occ, intercept, ...)
  barplot_inner(quantiles_p, names_detect, intercept, ...)
  
}

#' @export
#' @rdname methods
#' 
plot_pr_occ <- function(object, npred = 1000, var_name = NULL, label = NULL, ...) {
  
  if (is.null(var_name)) {
    var_idx <- 2
  } else {
    var_idx <- match(var_name, object$data$var_names_occ)
    var_idx <- var_idx + 1
  }
  
  if (is.null(label))
    label <- object$data$var_names_occ[var_idx]
  
  seq_occ <- array(0, dim = c(npred, dim(object$data$X_occ)[2]))
  seq_occ[, 1] <- seq_occ[, 1] + 1
  seq_detect <- array(0, dim = c(npred, dim(object$data$X_detect)[2:3]))
  seq_detect[, 1, ] <- seq_detect[, 1, ] + 1
  seq_occ[, var_idx] <- seq(min(object$data$X_occ[, var_idx], na.rm = TRUE),
                            max(object$data$X_occ[, var_idx], na.rm = TRUE),
                            length = npred)
  seq_data <- list(X_occ = seq_occ, X_detect = seq_detect)
  out <- predict(object, newdata = seq_data, type = "response")
  plot(out$p_occ ~ seq_occ[, var_idx],
       bty = "l", xlab = paste0("Standardised ", label), ylab = "Pr(occupancy)",
       las = 1, type = "l", lwd = 2, col = "black",
       ylim= c(0, 1))
  points(c(object$data$y) ~ rep(object$data$X_occ[, var_idx], 2), pch = 16,
         col = ggplot2::alpha("black", 0.5))
  
}

#' @export
#' @rdname methods
#' 
plot_pr_detect <- function(object, npred = 1000, var_name = NULL, label = NULL, scale = NULL, ...) {
  
  if (is.null(var_name)) {
    var_idx <- 2
  } else {
    var_idx <- match(var_name, object$data$var_names_detect)
    var_idx <- var_idx + 1
  }
  
  if (is.null(label))
    label <- object$data$var_names_detect[var_idx]
  
  if (is.null(scale)) {
    label <- paste0("Standardised ", label)
    scale <- c(0, 1)
  }
  
  seq_occ <- array(0, dim = c(npred, dim(object$data$X_occ)[2]))
  seq_occ[, 1] <- seq_occ[, 1] + 1
  seq_detect <- array(0, dim = c(npred, dim(object$data$X_detect)[2:3]))
  seq_detect[, var_idx, 1] <- seq(min(object$data$X_detect[, var_idx, 1], na.rm = TRUE),
                                  max(object$data$X_detect[, var_idx, 1], na.rm = TRUE),
                                  length = npred)
  seq_detect[, 1, ] <- seq_detect[, 1, ] + 1
  seq_data <- list(X_occ = seq_occ, X_detect = seq_detect)
  out <- predict(object, newdata = seq_data, type = "response")
  x_vals <- scale[1] + scale[2] * seq_detect[, var_idx, 1]
  plot(out$p_detect[, 1] ~ x_vals,
       bty = "l", xlab = label, ylab = "Pr(detection)",
       las = 1, type = "l", lwd = 2, col = "black",
       ylim= c(0, 1))
  x_pts <- scale[1] + scale[2] * object$data$X_detect[, var_idx, 1]
  points(c(object$data$y) ~ rep(x_pts, 2), pch = 16,
         col = ggplot2::alpha("black", 0.5))
  
}
