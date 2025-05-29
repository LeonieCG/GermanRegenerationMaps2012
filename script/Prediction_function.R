# These are the functions needed to predict species distribution maps
# NOTE: These functions are highly specific and only for bams with families using log link functions


# Create multi-layer raster with all predictors--------------------------------------------------------------
predict.stack <- function(vars,
                          exclude = NULL,
                          get_dir,
                          prefix = NULL,
                          aggreg = NULL){
    Start <- Sys.time()
    
    ifelse(is.null(exclude),
      vars <- vars,
      vars <- vars[!vars %in% exclude])
    
    vars <- case_match(vars,
                       "alt" ~ "alt_eudem", 
                       "northexp" ~ "northexp_der",
                       "eastexp"~ "eastexp_der", 
                       "blname" ~ "blname_bkg",
                       .default=vars) # variables that originate from original BWI data have an external raster source and therefore different name
    
    # Load predictor rasters of interest 
    Newstack <- rast(paste0(get_dir, prefix, vars[1], ".tif"))
    
    for (i in 2:length(vars)){
      stackadd <- rast(paste0(get_dir, prefix, vars[i], ".tif"))
      Newstack <- c(Newstack, stackadd)
      rm(stackadd)
    }
    
    # Aggregate
    ifelse(is.null(aggreg),
           Newstack <- Newstack,
           Newstack <- terra::aggregate(Newstack, fact = aggreg))# aggreg = 10 meaning x km resolution? default is building the mean
    
    # Rename columns
    Newstack %<>%
      tidyterra::rename(alt = alt_eudem,
                        northexp = northexp_der,
                        eastexp = eastexp_der) # see argumentation above
    
    return(Newstack)
  print(Sys.time()- Start)
}
  


# Prediction function for predicting on mulitlayer raster predictor data -----------------------------------------------------
predict.fit <- function(fit,
                        predictorstack,
                        exclude,
                        BAspecies = NULL, # otherwise give species (use for regeneration prediction, will add basal area for regarding species)
                        CI = FALSE, # if confidence intervals needed or not
                        haconvert = FALSE,
                        get_dir = NULL,
                        cores){
  Start <- Sys.time()

  # add wzp_ba_ha_species for the regarding tree species to the predictor raster dataset
  if(is.null(BAspecies)){
    Newstack <- predictorstack
    } else {
      ba <- rast(paste0(get_dir, "wzp12_ba_ha_species_", BAspecies, ".tif"))
      names(ba) <- "wzp12_ba_ha_species"
      Newstack <- c(predictorstack, ba)
      }

  # define important stuff again so it is in the functions environment (otherwise it does not work..)
  exclude <- exclude
  
  # Predict on parallel clusters
  cls <- parallel::makeCluster(cores) # set up cores for parallelization
  parallel::clusterExport(cls, c("fit", "predict.parallel", "predict.bam"), envir = environment()) # "exclude" needed when run directly 
  predict <- terra::predict(Newstack, fit, na.rm =TRUE, exclude = exclude, cls = cls, # terra.predict turns rasters into format that can be processed by predict.bam, in this case it is transferred to the function for predict.bam in parallel
                            fun = predict.parallel)
  parallel::stopCluster(cls) # stop clusters
  
  # Prediction is NOT on the response scale, use inverse of link function to have the original scale
  predict$expfit <- exp(predict$fit) # link is log therefore exp, inverse function is stored in fit$family$linkinv() can be used for inverse, but does not work in this case
  
  # Create confidence interval limits upper and lower if CI == TRUE and rename layers
  if(CI==TRUE){
    predict$upperCI <- exp(predict$fit + (1.96 * predict$se.fit))
    predict$lowerCI <- exp(predict$fit - (1.96 * predict$se.fit))

    predict %<>%
      select(c(expfit,upperCI,lowerCI)) %>%
      rename(!!attr(fit, "species") := expfit)
  } else {
    predict %<>%select(c(expfit)) %>%
      rename(!!attr(fit, "species") := expfit)
  }

  # Conversion from count per 2 radius sample plots to count per ha
  ifelse(haconvert == FALSE,
         predict <- predict,
         predict <- (predict*10000)/(2^2*pi))
  
  # Return result
  return(predict)

  print(Sys.time()- Start)
}


# Predict in parallel --------------------------------------------------------
predict.parallel <- function(mod, dat, cls, exclude, ...) {
  ncls <- length(cls) # number of cores
  nr <- nrow(dat)# dataset rows
  splitted <- split(dat, rep(1:ncls, each=ceiling(nr/ncls), length.out=nr))  # split dataset in chunks for each core
  
  #take one chunk for each core and predict
  predictlist = parallel::clusterApply(cls, splitted,
                                       function(x,...) predict.bam(mod, x,
                                                                   # type= "response", # for predicting on response scale, if turned off predict needs exp(predict) due to the link function
                                                                   exclude = paste0("s(", exclude, ")"),
                                                                   newdata.guaranteed = TRUE,
                                                                   discrete = FALSE,
                                                                   se.fit = TRUE, # for return of standard errors!
                                                                   ...))
  
  #each core returns $fit and $se.fit and stores it in list per core, $fit and $se.fit need to be combined across all cores and passed back to terra::predict
  list(fit=unname(unlist(lapply(predictlist, function(x) c(x$fit)))),
       se.fit=unname(unlist(lapply(predictlist, function(x) c(x$se.fit)))))# returns list with $fit and $fit.se
}

