# internal wrapper for mean function to handle all NA obs
mean_fn <- function(x, ...) {
  
  out <- NA
  if (!all(is.na(x)))
    out <- mean(x, ...)
  
  out
  
}
