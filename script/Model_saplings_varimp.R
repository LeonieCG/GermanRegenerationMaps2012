# This script calibrates species-specific regeneration models while leaving out variable groups
# this variable group importance is tested only for final.species

start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", "DHARMa",
              "mgcViz", "tidyterra", "parallel", "cowplot",
              "modelr", "sf","blockCV", "automap", "stringr",
              "parallel")
# sapply(packages, FUN = install.packages, character.only = T)
sapply(packages, FUN = library, character.only = T)


# System setup for parallel -----------------------------------------------
Sys.setenv(OMP_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           RHPCBLAS_CTL = "1")

# Choose how many cores to use for species-level parallelism
mc_cores_default <- max(1, parallel::detectCores() - 1)


# Load data ---------------------------------------------------------------

## Response Variable -------------------------------------------------------
regeneration <- readRDS("data/DE_BWI3_regeneration_h50d7.rds")

## Explanatory variables ---------------------------------------------------
# Species independent
E <- readRDS("data/DE_BWI3_explvars_df.rds") 

# Species dependent
WZP12 <- readRDS("data/DE_BWI3_big_basalarea_wzp12_interpol.rds")

## Sapling species -----------------------------------------------------------------
species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species

## Source functions --------------------------------------------------------
source("script/Model_functions.R")  # includes building Model variable table, Model function, Cross validation

## Choose variables ----------------------------------------------------------
# # due to licencing, not all environmental predictors could be published with this study 
# # the following is a selection, model outcomes can thereby vary highly from the study!
# resp = "count"
# fixed = c("tPeriodic2020_forcli", "tSeas2020_forcli","tMinColdMonth2020_forcli","tRangeDay2020_forcli","tRangeAn2020_forcli" #T microclimate
#           , "tPeriodic2010_chelsa", "tSeas2010_chelsa","tMinColdMonth2010_chelsa","tRangeDay2010_chelsa","tRangeAn2010_chelsa" #T macroclimate
#           , "precPeriodic2010_chelsa", "precSeas2010_chelsa" # prec
#           , "wwpi_cop" # Water prob index, Anoxy indicator, water bodies and flooded plains
#           , "tcd_cop" #tree cover density
#           , "alt", "northexp", "eastexp" # terrain vars
#           , "wzp12_ba_ha_species" #Basal area of respective old trees
# )
# random = c("yearmonth" # to account for changes within sampling period
#            , "blname" # Bundesland
# )
# spatial = c("x","y")

# Variables of the full model
source("script/Model_vars.R")

resp.full = resp
fixed.full = fixed
random.full = random
spatial.full = spatial 

rm(resp, fixed, random, spatial) # to avoid any mistakes


## Specify explanatory and prediction variables
resp = "count"
random = c("yearmonth" # year and month of measuring period
           , "blname") # Bundesland)
spatial = c("x","y") # xy coordinates of inspire-grid

varimp <- read.csv2("data/Model_vars_varimp.csv") %>% 
  filter(!(Category %in% c("Space", "Time")))# leave out these categories!

# categories that are leftout
leftouts <- unique(varimp$Category)


# ---- Main loop over leftout (serial) -------------------------------------
for (leftout in leftouts) {
  message(sprintf("left out: %s", leftout))
  
  # define fixed effects excluding current category
  fixed <- varimp %>%
    dplyr::filter(!(Category == leftout)) %>%
    dplyr::pull(varid)
  
  # create folder structure for this leftout (pre-create to avoid races)
  out_root <- file.path("output", "Fits", "Sapling", paste0("h50d7_Germany_varimp_", leftout))
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  for (species in species.final) {
    dir.create(file.path(out_root, species), recursive = TRUE, showWarnings = FALSE)
  }
  
  # Decide cores per leftout/species batch (cap by number of species to save RAM)
  mc_cores <- min(length(species.final), mc_cores_default)
  
  # Reproducibility across forks (important if CV/randomness is used)
  set.seed(84)
  
  # ---- Parallelize over species ------------------------------------------
  results <- mclapply(
    X = species.final,
    FUN = function(species) {
      # Worker-side progress
      message(sprintf("[pid=%s] species=%s | leftout=%s", Sys.getpid(), species, leftout))
      
      # Build model variables
      mv <- try({
        mv0 <- modelVariables(species = species)
        mv0 %>%
          dplyr::select(dplyr::all_of(c(
            resp.full, fixed.full, random.full, spatial.full,
            "plotid", "clusterid"
          ))) %>%
          tidyr::drop_na()
      }, silent = TRUE)
      
      if (inherits(mv, "try-error")) {
        warning(sprintf("modelVariables/select failed for species=%s; skipping.", species))
        return(list(species = species, status = "modelVariables_failed"))
      }
      
      # Fit the model (robust error handling)
      fit <- try({
        model.fit(resp = resp,
                  fixed = fixed,
                  fixedfact = NULL,
                  random = random,
                  spatial = spatial,
                  offset = NULL,
                  exclude = "yearmonth",
                  fam = nb,
                  s.k = 10,
                  s.bs = "cs",
                  spat.which = "tensor",
                  spat.bs = "cs", # is faster than ts
                  spat.k= "c(25,50)", 
                  select.var = FALSE,
                  bam = TRUE,
                  Data = mv,
                  CV = "blockcv",
                  blockcv.dr = 300000, # 300000 gives 11 blocks for whole germany
                  blockcv.k = 10
        )
      }, silent = TRUE)
      
      if (inherits(fit, "try-error")) {
        warning(sprintf("model.fit failed for leftout=%s species=%s; skipping save.", leftout, species))
        return(list(species = species, status = "model_fit_failed"))
      }
      
      # Save result
      save_path <- file.path(out_root, species, paste0(species, "_fit.rds"))
      try(saveRDS(fit, save_path), silent = TRUE)
      
      # Free memory in the worker
      rm(fit, mv); gc()
      
      list(species = species, status = "ok")
    },
    mc.cores = mc_cores,
    mc.preschedule = FALSE,  # better load balancing when species vary in duration
    mc.set.seed = TRUE       # independent RNG streams derived from set.seed above
  )
  
  # Optional: inspect per-leftout results summary
  ok_count <- sum(vapply(results, function(x) identical(x$status, "ok"), logical(1)))
  message(sprintf("leftout=%s done. %d/%d species completed.",
                  leftout, ok_count, length(species.final)))
  
  # Model checks-------------------------------------------------------------
  
  ## Start Summary over all species ------------------------------------------------
  # Save Model Summary in all species data frame
  df.allsp <- data.frame()
  
  ## Iterate Species ---------------------------------------------------------
  for(species in species.final) {
    print(species)
    
    try({
      # set up output chart
      df.out <-  data.frame()
      
      #Load fit
      fit <- readRDS(paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_fit.rds"))
      
      # Save Model output in data.frame
      df.out <- bind_rows(df.out, fit$CV[c("cv.method", "cv.folds", "cv.sp.range", "cv.set.range", "cv.blocknr")])
      
      
      ## Model validation indicator ----------------------------------------------
      # MAE
      df.out$mae.train.mean = mean(fit$CV$mae.train)
      df.out$mae.train.median = median(fit$CV$mae.train)
      df.out$mae.train.sd = sd(fit$CV$mae.train)
      df.out$mae.train.iqr = quantile(fit$CV$mae.train, 0.75) - quantile(fit$CV$mae.train, 0.25)
      
      df.out$mae.test.mean = mean(fit$CV$mae.test)
      df.out$mae.test.median = median(fit$CV$mae.test)
      df.out$mae.test.sd = sd(fit$CV$mae.test)
      df.out$mae.test.iqr = quantile(fit$CV$mae.test, 0.75) - quantile(fit$CV$mae.test, 0.25)
      
      df.out$mae.relative.mean = mean(fit$CV$mae.test/fit$CV$mae.train)
      df.out$mae.relative.median = median(fit$CV$mae.test/fit$CV$mae.train)
      
      #Cohens pseudo R2
      df.out$rsq.train.mean = mean(fit$CV$rsq.train)
      df.out$rsq.train.median = median(fit$CV$rsq.train)
      df.out$rsq.train.sd = sd(fit$CV$rsq.train)
      df.out$rsq.train.iqr = quantile(fit$CV$rsq.train, 0.75) - quantile(fit$CV$rsq.train, 0.25)
      
      df.out$rsq.test.mean = mean(fit$CV$rsq.test)
      df.out$rsq.test.median = median(fit$CV$rsq.test)
      df.out$rsq.test.sd = sd(fit$CV$rsq.test)
      df.out$rsq.test.iqr = quantile(fit$CV$rsq.test, 0.75) - quantile(fit$CV$rsq.test, 0.25)
    })# close try
    
    ## Data stats --------------------------------------------------------------
    Data <- modelVariables(species = species) %>%
      select(all_of(c(resp, fixed, random, spatial))) %>%
      drop_na()
    
    # How many plots are left when Nas dropped?
    df.out[1, "data.n"] <- dim(Data)[1]
    
    # How many plots are containing zeros?
    df.out[1, "data.resp.n0"] <- Data %>% select(all_of(resp)) %>% summarise(sum(.==0))
    
    
    ## Save outputs -------------------------------------------------------------
    #df.out
    df.out[1,"species"] <- species
    write.csv(df.out, paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_fit.csv"), row.names = F)
    
    # save in all species df
    df.allsp <- bind_rows(df.allsp, df.out[1,])
  }
  
  
  ## Save summary over all species -------------------------------------------
  # Save df
  write.csv2(df.allsp, paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/Sapling_model_summary.csv"), row.names = FALSE)
}



print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))