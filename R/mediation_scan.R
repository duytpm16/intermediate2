#'
#' Mediation Scan
#' Scan set of mediators for a target.
#' 
#' @param target A numeric vector with gene/protein expression
#' @param mediator A matrix, each column is one gene/protein's expression
#' @param driver A matrix, haplotype probabilities at QTL we try to mediate
#' @param annotation A data frame with mediators' annotation with columns `facet_name` and `index_name`
#' @param covar A matrix with additive covariates
#' @param intcovar A matrix of covariate interacting with driver
#' @param kinship kinship object
#' @param method A method to handle missing cases
#' @param fitFunction function to fit models
#' @param facet_name name of facet column (default `chr`)
#' @param index_name name of index column (default `pos`)
#' @param verbose If TRUE display information about the progress
#' 
#' @details 
#' For a given QTL haplotype probabilities `driver`` and target `target`,
#' the function sequentially tries to add each column of `mediator` matrix as a covariate
#' and calculates LOD statistic. The low LOD value indicates `driver` and
#' `target` are conditionally independent given `mediator`,
#' i.e. `mediator` is a mediator of causal relationship from `driver` to `target`.
#'
#' @examples
#' data(Tmem68)
#' target <- Tmem68$target
#' m <- match("Tmem68", Tmem68$annotation$symbol)
#' 
#' med_scan <- mediation_scan(target = target,
#'                       mediator = Tmem68$mediator,
#'                       driver = Tmem68$qtl.geno,
#'                       annotation = Tmem68$annotation,
#'                       covar = Tmem68$covar,
#'                       method = "double-lod-diff")
#' ggplot2::autoplot(med_scan)
#' ggplot2::autoplot(subset(med_scan, "4")) +
#'   ggplot2::geom_vline(xintercept = Tmem68$annotation[m,"pos"], linetype = "dashed")
#' 
#' @export

mediation_scan <- function(target, mediator,  driver, 
                           annotation, covar=NULL, intcovar=NULL,
                           kinship = NULL, method=c("double-lod-diff", "ignore", "lod-diff"), 
                           fitFunction = fitDefault, facet_name = "chr", index_name = "pos", verbose=TRUE) {
  
  
  

        # Make sure covariates are numeric
        covar    <- covar_df_mx(covar)
        intcovar <- covar_df_mx(intcovar)
        
        
        
        
        
        
        # Fit model without mediation.
        if(length(dim(driver)) == 2){
           driver_tar <- driver
        }else{
           driver_tar <- driver[,,1]
        }

        loglik0 <- fitFunction(driver_tar, target, kinship, covar, intcovar)$LR
        
        
      
  
        # Check input
        stopifnot(c(facet_name, index_name) %in% tolower(names(annotation)))
  
        # Match up annotation with mediators
        stopifnot(all(!is.na(m <- match(colnames(mediator), annotation$id))))
        annotation <- annotation[m,]
        rownames(annotation) <- colnames(mediator)
        

        
  
  
        
        mapfn <- function(x, target, covar, driver, loglik0) {
        
          
                # Get common data.
                commons <- common_data(target = target, mediator = x$mediator, driver = driver, covar = covar, intcovar = intcovar)
        
                if(is.null(commons)){
                   return(NULL)
                }
        
                
                target   <- commons$target
                mediator <- commons$mediator
                driver   <- commons$driver
                covar    <- commons$covar_tar
                intcovar <- commons$intcovar
                common   <- commons$common
                rm(commons)
        
                
                
                method = match.arg(method)
        
      
                lodfn <- switch(method,
                                ignore            = function(loglik, loglik0) loglik[1],
                                "lod-diff"        = function(loglik, loglik0) loglik[2] - loglik[1],
                                "double-lod-diff" = function(loglik, loglik0) loglik0 - (loglik[2] - loglik[1]))
        
      
                
                
                
                if(length(dim(driver)) > 2){
                   if(is.null(dcol <- x$annotation$driver_names))
                      driver <- driver[,,1]
                   else
                      driver <- driver[,, dcol]
                }
        
        
                loglik <- c(fitFunction(driver, target, kinship,
                                        cbind(covar, mediator), intcovar)$LR,
                            fitFunction(driver, target, kinship,
                                        covar, intcovar)$LR)
                
                
                lodfn(loglik, loglik0)
        } # mapfn()
        
        
                                
                                
                                
        # List of mediators
        med_pur <- purrr::transpose(list(mediator = as.data.frame(mediator),
                                         annotation = split(annotation, rownames(annotation))))
        
        
        # Add mediation LOD to annotation
        output     <- annotation
        output$lod <- unlist(purrr::map(med_pur, mapfn, target, covar, driver, loglik0)) / log(10)
        
        
        
        attr(output, "targetFit")  <- loglik0 / log(10)
        attr(output, "facet_name") <- facet_name
        attr(output, "index_name") <- index_name
        class(output) <- c("mediation_scan", "data.frame")
        return(output)
}






#'
#' @export
subset.mediation_scan <- function(object, facets=NULL, ...) {
        facet_name <- attr(object, "facet_name")
        
        
        if(is.null(facets)){
           return(object)
        }
        
        
        new_object <- object[object[[facet_name]] %in% facets,]
        modify_object(object, new_object)
}
