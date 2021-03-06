#' @rdname fitDefault
#' @export
fitQtl2 <- function(driver,
                    target,
                    ...) {

  if(is.null(rownames(driver)))
    rownames(driver) <- seq_len(nrow(driver))
  
  out <- qtl2::fit1(driver, target, ...)
  
  # Replace lod names with LR
  names(out) <- stringr::str_replace(names(out), "lod", "LR")
  names(out) <- stringr::str_replace(names(out), "_LR", "LR")

  # Rescael to make them likelihoods (or likelihood ratios)
  out$LR <- out$LR * log(10)
  out$indLR <- out$indLR * log(10)
  
  # Add df for later use
  out$df <- ncol(driver) - 1
  
  # Residuals
  fitted <- rep(NA, length(target))
  names(fitted) <- if(is.matrix(target)) {
    rownames(target)
  } else {
    names(target)
  }
  fitted[names(out$fitted)] <- out$fitted
  out$resid <- target - fitted
  
  out
}