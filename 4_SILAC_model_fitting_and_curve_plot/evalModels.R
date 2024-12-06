evalModels <- function(fits, groups, x){
  # Invoke model evaluation
  
  modelMetrics <- assessModelMetrics(fits = fits, x = x, groups = groups)
  modelPredictions <- augmentModels(fits = fits, groups = groups)
  
  return(list(modelMetrics = modelMetrics, modelPredictions = modelPredictions))
}

augmentModels <- function(fits, groups){
  # Augment model to obtain data frame with measurements, predictions and residuals.
  
  # ---- Computing model predictions and residuals ----
  message("Computing model predictions and residuals ...")
  
  modelPredictions <- fits %>%
    filter(.data$successfulFit) %>%
    group_by_at(c("iter", "id", groups)) %>%
    do(augment(.data$fittedModel[[1]])) %>%
    ungroup %>%
    right_join(fits %>% select(!!c("iter", "id", groups)), by = c("iter", "id", groups)) %>% # fill for unsuccessful fits
    select(-.data$iter)
  
  message("... complete.\n")
  
  return(modelPredictions)
  
}

assessModelMetrics <- function(fits, x, groups){
  # Invoke calculation of performance metrics per model
  
  message("Evaluating models ...")
  
  metrics <- fits %>%
    group_by_at(c("iter", "id", groups)) %>%
    do(assessSingleModel(nls_obj = .data$fittedModel[[1]],
                         xVec = unique(x))) %>%
    ungroup %>%
    select(-.data$iter)
  
  message("... complete.\n")
  
  return(metrics)
}

#' @importFrom methods is
assessSingleModel <- function(nls_obj, xVec){
  
  if (any(is.na(xVec))){
    stop("Temperature vector contains missing values. Cannot integrate over these values to compute the curve area.")
  }
  
  tm = a = b = aumc = tm_sd = rss = logL = nCoeffs = nFitted = resid_sd <- NA
  conv <- FALSE
  if (!is(nls_obj, "try-error")){
    #meltCurve <- Reduce(paste, deparse(formula(nls_obj))) %>% gsub("y ~ ", "", .)
    pars <- coefficients(nls_obj)
    nCoeffs <- length(pars)
    a <- pars[["a"]]
    b <- pars[["b"]]
    yVec <- fitted(nls_obj) + resid(nls_obj)
    tm_sd <- summary(nls_obj)$coefficients[, "Std. Error"][2]
    ## Compute area under the melting curve:
    int <- try(integrate(function(x,m) predict(m, newdata = list(x=x)),
                         min(xVec), max(xVec), nls_obj), silent = TRUE)
    if (!inherits(int, "try-error")){
      aumc <- int$value
    }
    
    nFitted <- nobs(nls_obj)
    rss <- sum(resid(nls_obj)^2, na.rm = TRUE)
    resid_sd  <- sqrt(rss/nFitted)
    logL <- -nFitted/2 * log(2*pi*resid_sd^2) - rss/(2*resid_sd^2) #loglik <- logLik(m)
    nonNAs <- sum(!is.na(resid(nls_obj)))
    
    conv <- nls_obj$convInfo$isConv
  }
  return(data.frame(tm = b, a = a,  b = b, aumc = aumc,
                    resid_sd = resid_sd, rss = rss, loglik = logL,
                    tm_sd = tm_sd, nCoeffs = nCoeffs, nFitted = nFitted, conv = conv))
}