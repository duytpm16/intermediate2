common_data <- function(target = NULL, mediator = NULL, driver = NULL,
                        covar_tar = NULL, covar_med = NULL, kinship = NULL,
                        driver_med = NULL, intcovar = NULL,
                        common = TRUE,
                        minN = 0, minCommon = 0.9, ...) {

  # Make sure all are matrices
  target <- convert_matrix(target, "T")
  driver <- convert_matrix(driver, "D", rownames(target))
  mediator <- convert_matrix(mediator, "M", rownames(target))
  
  covar_tar <- convert_matrix(covar_tar,
                              paste0("covT", seq_len(ncol(covar_tar))), 
                              rownames(target))
  covar_med <- convert_matrix(covar_med,
                              paste0("covM", seq_len(ncol(covar_med))), 
                              rownames(mediator))
  driver_med <- convert_matrix(driver_med,
                              paste0("driverM", seq_len(ncol(driver_med))), 
                              rownames(mediator))
  intcovar <- convert_matrix(intcovar,
                             paste0("intcov", seq_len(ncol(intcovar))), 
                             rownames(target))
  
  if(is.list(kinship)) {
    if(is.matrix(kinship[[1]]))
      kinship <- kinship[[1]]
  }
  K <- kinship  
  if(!is.matrix(K))
    K <- NULL
  
  # Keep individuals with full records.
  ind2keep <- get_common_ids(driver, target, covar_tar, covar_med, K, driver_med, intcovar,
                             complete.cases = TRUE)
  
  # Drop mediator columns with too few non-missing data.
  if(enough <- (length(ind2keep) >= 0)) {
    m <- match(ind2keep, rownames(mediator), nomatch = 0)
    ind2keep <- ind2keep[m > 0]
    mediator <- mediator[m,, drop = FALSE]

    # This way considers only ind with no missing data.
    # Might want another way if pattern of missing different for some mediators.
    # Then would not do decomp_kinship and need to do common_data for single mediator
    
    # Count as number if not infinite and not missing
    is_num <- function(x) { !is.na(x) & is.finite(x) }
    # Drop mediators with too little data.
    ok_med <- apply(mediator, 2, function(x) sum(is_num(x))) >= 0
    if(enough <- any(ok_med)) {
      mediator <- mediator[, ok_med, drop = FALSE]

      # Check for missing across all remaining mediators.
      allMed <- apply(mediator, 1, function(x) any(is_num(x)))
      mediator <- mediator[allMed,, drop = FALSE]
      ind2keep <- ind2keep[allMed]
      if(enough <- (length(ind2keep) >= 0)) {
        common <- common & (sum(ok_med) == 1) | all(is_num(mediator))

      }
    }
  }
  if(!enough) {
    warning(paste0("too few data (", length(ind2keep), ") for analysis"))
    return(NULL)
  }
  driver <- switch(
    length(dim(driver)),
    as.matrix(driver)[ind2keep,, drop = FALSE],
    driver[ind2keep,, drop = FALSE],
    driver[ind2keep,,, drop = FALSE]
  )
  target <- target[ind2keep,, drop = FALSE]
  mediator <- mediator[ind2keep,, drop = FALSE]
  if(!is.null(covar_tar))
    covar_tar <- covar_tar[ind2keep,, drop = FALSE]
  if(!is.null(covar_med))
    covar_med <- covar_med[ind2keep,, drop = FALSE]
  if(!is.null(driver_med))
    driver_med <- driver_med[ind2keep,,, drop = FALSE]
  if(!is.null(intcovar))
    intcovar <- intcovar[ind2keep,, drop = FALSE]
  
  if(!is.null(kinship)) {
    if(is.matrix(kinship))
      kinship <- kinship[ind2keep, ind2keep]
    # Decompose kinship if all in common.
    is_kinship <- attr(kinship, "eigen_decomp")
    if(is.null(is_kinship))
      is_kinship <- FALSE
    if(common && !is_kinship)
      kinship <- qtl2::decomp_kinship(kinship)
  }
  list(driver = driver,
       target = target,
       mediator = mediator,
       kinship = kinship,
       covar_tar = covar_tar,
       covar_med = covar_med,
       driver_med = driver_med,
       intcovar = intcovar,
       common = common)
}
convert_matrix <- function(object, 
                           col_names = "",
                           row_names = seq_len(nrow(object))) {
  if(is.null(object))
    return(NULL)
  
  if(!is.array(object))
    object <- as.matrix(object)
  stopifnot(is.numeric(object))
  if(is.null(colnames(object)) & !is.null(col_names))
    colnames(object) <- col_names
  if(is.null(rownames(object)) & !is.null(row_names))
    rownames(object) <- row_names
  object
}
