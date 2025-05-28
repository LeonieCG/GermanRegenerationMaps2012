# Model Variables ---------------------------------------------------------------
modelVariables <- function(species){
  R <- as.data.frame(regeneration$data[regeneration$tax == species])
  B12 <- as.data.frame(WZP12$data[WZP12$tax == species])# Basal area of wzp12
  
  # join Saplings, Basal Area and Environmental Variables
  mv <- R %>% 
    dplyr::left_join(B12, by = "plotid") %>% 
    dplyr::left_join(E,., by = "plotid")

  mv %<>% mutate(clusterid = as.factor(clusterid))
  
  # add year month variable
  mv %<>%
    mutate(yearmonth =  as.factor(format(as.Date(mv$time,format ="%Y-%m-%d"),"%Y-%m")))

  mv %<>% rename(alt = alt_loc,
                 northexp = northexp_loc,
                 eastexp = eastexp_loc,
                 blname = blname_loc,
                 ownership = ownership_loc)  
  
  attr(mv, "species") <- species # get information via attributes(mv)$species
  return(mv)
}



# Model function ----------------------------------------------------------

model.fit <- function(resp, 
                      fixed = NULL,
                      fixedfact = NULL,
                      random = NULL,
                      spatial = NULL,
                      offset = NULL,
                      exclude = NULL,
                      s.k = -1 , # Default s()
                      te.k = NA, # Default te()
                      ste.bs = "tp",
                      select.var = FALSE,
                      fam,
                      bam = FALSE,
                      Data,
                      CV = FALSE, # give blockcv or False = no cv
                      blockcv.dr = NULL,# blockcv default range
                      blockcv.k # Number of folds
) {
  start.model <- Sys.time()
  
  coordinates <- c("x","y") # add x,y coordinates also if spatial = NULL because it is needed for spat. auto corr range calc of blocked cv
  
  # select needed columns
  Data %<>% 
    select(all_of(c(resp, fixed, fixedfact, random, coordinates))) %>% 
    drop_na()
  
  # built formula parts and whole formula
  form_resp = ifelse(length(resp) >=2, paste0("cbind(", paste(resp, collapse = ","), ")"), paste0(resp))
  form_fixed <-  ifelse(is.null(fixed) , paste0("1"), paste(paste0("s(", fixed, ", bs= '", ste.bs, "', k = ", s.k,")"), collapse = " + "))
  form_fixedfact  <- ifelse(is.null(fixedfact) , paste0(""), paste(paste0("+",fixedfact), collapse = " "))
  form_random <- ifelse(is.null(random) , paste0(""), paste(paste0("+ s(",random, ", bs='re')"), collapse = " "))
  form_spatial <- ifelse(is.null(spatial) , paste0(""), paste0("+ te(", paste0(spatial, collapse = ","), ", bs= '", ste.bs, "',  k = ", te.k,")"))
  form_offset <- ifelse(is.null(offset) , paste0(""), paste(paste0("+ offset(log(",offset, "))"), collapse = " "))
  
  form_all <- as.formula(paste(form_resp," ~ ",form_fixed, form_fixedfact, form_random, form_spatial, form_offset))
  print(form_all)
  
  # Bam or gam
  if(bam == FALSE){#gam
    fit <- gam(formula = form_all, family = fam, data = Data
               , select = select.var
               , method = "REML") # method: method to choose lambda (weighing of wiggliness), setting the penalty for complexity (Gavin Simpson video)
  } else {#bam
    fit <- bam(formula = form_all, family = fam, data = Data
               , select = select.var
               , method = "fREML" 
               , discrete=TRUE, nthreads=10) # speeds up calculation
  }
  print(paste("Model took", Sys.time() - start.model, units(Sys.time()-start.model)))
  
  # Cross validation
  if(CV == "blockcv"){
    fit$CV <- blockcv(blockcv.k = blockcv.k,
                      blockcv.dr = blockcv.dr,
                      spatial = spatial,
                      select.var = select.var,
                      Data = Data,
                      form_all = form_all,
                      fam = fam,
                      exclude = exclude,
                      resp = resp,
                      bam = bam)
  } else {
    print("CV not applied")
  }
  
  attr(fit, "species") <- attr(Data, "species")
  print(paste("Model + CV took", Sys.time() - start.model, units(Sys.time()-start.model)))
  
  return(fit)
}




# Blocked cross validation ------------------------------------------------

packages <- c("modelr", "sf","blockCV", "automap") 
sapply(packages, FUN = library, character.only = T)

blockcv <- function(blockcv.k,
                    blockcv.dr = NULL,
                    Data,
                    form_all,
                    spatial,
                    select.var,
                    fam,
                    exclude,
                    resp,
                    bam){               
  start.cv <- Sys.time()
  
  # Data
  Data_sf <-  Data %>% 
    sf::st_as_sf(., coords = c("x", "y"), crs = crs(readRDS("data/DE_BWI3_explvars_sv.rds")))
  
  # Check spatial autocorrelation distance in raw data
  range <- cv_spatial_autocor(x = Data_sf, # data 
                              column = resp, # response column
                              plot = FALSE,
                              progress = FALSE) # turns of folds
  
  # Choose distance used to built blocks:
  # if null use spatial autocorr range
  # if species spatial autocorrelation range bigger than threshold, set threshold
  # because then blocks <10, but min 10 needed for 10-fold block cv
  if(is.null(blockcv.dr)){
    cv.range <- ceiling(range$range)
  } else if(ceiling(range$range) <= blockcv.dr){
    cv.range <- ceiling(range$range)
  } else {
    cv.range <- blockcv.dr
  }
  
  # Set up test and training data subsets
  set.seed(15) # Set seed for comparability of results
  scv <- cv_spatial(x = Data_sf,
                    k = blockcv.k, # number of folds
                    size = cv.range, # size of the blocks in meters selected via autocorrelation function
                    selection = "random", # random blocks-to-fold
                    plot = FALSE,
                    progress = FALSE) 
  
  
  # perform cross validation iterating through all folds 
  mae.train <- vector()
  mae.test <- vector()
  rsq.train <- vector()
  rsq.test <-  vector()
  
  # if spatial = NULL omit coordinates from Data
  if(is.null(spatial)){
    Data <- Data %>% select(-c(x,y))
  } # else coordinates remain
  
  for (i in 1:blockcv.k){
    # defining test and training data set
    train <- Data[scv$folds_list[[i]][[1]],]
    test <-  Data[scv$folds_list[[i]][[2]],]
    
    if(bam==FALSE){
      fit.cv <- gam(formula = form_all, family = fam, data = train,
                    method = "REML",
                    select = select.var)
      fit.null <- gam(formula = as.formula(paste(resp," ~ 1")), family = fam, data = train,
                      method = "REML",
                      select = select.var)
      prediction.cv <- predict.gam(fit.cv,
                                   newdata = test[,!names(test) %in% exclude],
                                   type = "response", 
                                   exclude = paste0("s(", exclude, ")"),
                                   newdata.guaranteed = TRUE, 
                                   na.rm = TRUE)
    } else {
      fit.cv <- bam(formula = form_all, family = fam, data = train,
                     method = "fREML", 
                     select = select.var,
                     discrete = TRUE, 
                    nthreads = 10)
      fit.null <- gam(formula = as.formula(paste(resp," ~ 1")), family = fam, data = train,
                      method = "REML",
                      select = select.var) # gam here since bam does not work for null model
      prediction.cv <- predict.bam(fit.cv, 
                                   newdata = test[,!names(test) %in% exclude],# don't provide data that is excluded
                                   type = "response", 
                                   exclude = paste0("s(", exclude, ")"),
                                   newdata.guaranteed = TRUE, #due to exclude
                                   na.rm = TRUE,
                                   discrete = FALSE)# predict.bam specific
    } 
    
    # MAE
    mae.train[i] <- mean(abs(train[, resp] - fit.cv$fitted.values), na.rm = T)
    mae.test[i] <- mean(abs(test[, resp] - prediction.cv), na.rm = T)
    
    #R-Squared
    rsq = get_cohenrsq(fit = fit.cv, fit_null = fit.null, train = train, test = test, resp = resp,  exclude = exclude)
    
    rsq.train[i] <- rsq$cohenrsq.train
    rsq.test[i] <- rsq$cohenrsq.test
    }
  
  
  # Save CV output with model information
  CV <- list(cv.method = "blocked cross validation", 
             cv.folds = blockcv.k, 
             cv.sp.range = range$range,
             cv.set.range = cv.range,
             cv.blocknr = dim(scv$blocks)[1],
             mae.train = mae.train,
             mae.test = mae.test,
             rsq.train = rsq.train,
             rsq.test = rsq.test)

  print(paste("CV took", Sys.time() - start.cv, units(Sys.time()-start.cv)))
  structure(CV)
}


# Pseudo R2 by Cohen ------------------------------------------------------

get_cohenrsq <- function(fit, 
                         fit_null, 
                         train, 
                         test, 
                         resp, 
                         exclude){
  
  model_type <- class(fit)[1] # determine model type
  family <- switch(model_type,
                   gam = strsplit(fit$family$family, '[()]|\\ ')[[1]][1],
                   bam = strsplit(fit$family$family, '[()]|\\ ')[[1]][1]) # retrieve distribution family
  
  if(is.null(family)|!family %in% c("Negative")) stop(paste(family, "family not supported"))
  if(!model_type %in% "bam") stop(paste(modeltype, "class not supported"))
  
  
  # Get mean
  y_mean <- mean(predict.bam(fit_null, 
                            newdata = train[,!names(train) %in% exclude],# don't provide data that is excluded
                            type = "response", 
                            exclude = paste0("s(", exclude, ")"),
                            newdata.guaranteed = TRUE, #due to exclude
                            na.rm = TRUE,
                            discrete = FALSE)) # predict.bam specific
  
  # Predicted response
  y_pred.train <- predict.bam(fit, 
                             newdata = train[,!names(train) %in% exclude],# don't provide data that is excluded
                             type = "response", 
                             exclude = paste0("s(", exclude, ")"),
                             newdata.guaranteed = TRUE, #due to exclude
                             na.rm = TRUE,
                             discrete = FALSE)
  
  y_pred.test <- predict.bam(fit, 
                            newdata = test[,!names(test) %in% exclude],# don't provide data that is excluded
                            type = "response", 
                            exclude = paste0("s(", exclude, ")"),
                            newdata.guaranteed = TRUE, #due to exclude
                            na.rm = TRUE,
                            discrete = FALSE)
  
  # Log-Likelihood
  LL_null.train = sum(dnbinom(train[,resp], size = fit_null$family$getTheta(TRUE), mu = y_mean, log = T)) # null
  LL_sat.train = sum(dnbinom(train[,resp], size = fit$family$getTheta(TRUE), mu = train[,resp], log = T)) # saturated
  LL_pred.train = sum(dnbinom(train[,resp], size = fit$family$getTheta(TRUE), mu = y_pred.train, log = T)) # pred
  
  LL_null.test = sum(dnbinom(test[,resp], size = fit_null$family$getTheta(TRUE), mu = y_mean, log = T)) # null
  LL_sat.test = sum(dnbinom(test[,resp], size = fit$family$getTheta(TRUE), mu = test[,resp], log = T)) # saturated
  LL_pred.test = sum(dnbinom(test[,resp], size = fit$family$getTheta(TRUE), mu = y_pred.test, log = T)) # pred
  
  # Deviance
  dev_pred.train = 2 * (LL_sat.train - LL_pred.train)
  dev_null.train = 2 * (LL_sat.train - LL_null.train)
  
  dev_pred.test = 2 * (LL_sat.test - LL_pred.test)
  dev_null.test = 2 * (LL_sat.test - LL_null.test)
  
  # Cohens pseudo R square
  cohenrsq.train = 1 - dev_pred.train/dev_null.train
  cohenrsq.test = 1 - dev_pred.test/dev_null.test
  
  # Return values
  structure(
    list(cohenrsq.train = cohenrsq.train,
         cohenrsq.test = cohenrsq.test)
  )
}

