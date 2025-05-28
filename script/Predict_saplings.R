start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", 
              "mgcViz", "tidyterra", "parallel")
sapply(packages, FUN = library, character.only = T)


# Choose location --------------------------------------------------
outloc = "h50d7_Germany/" # in this case location of fit


# Load data ---------------------------------------------------------------
# Sapling species
species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species


# Source functions --------------------------------------------------------
source("script/Prediction_function.R")


# Choose variables ----------------------------------------------------------
source("script/Model_vars.R")
fixed = fixed[! fixed %in% c("wzp12_ba_ha_species")] # removed here because it it will be added later predict.fit(BAspecies=species)


# Predictor stack ---------------------------------------------------------
stacked <-  predict.stack(vars =  c(fixed, spatial,random),
                          exclude = "yearmonth",
                          get_dir = "data/Predictor_100m_Germany/")


# Run prediction ----------------------------------------------------------
for(species in species.final) {
  print(species)
  try({
    prediction = predict.fit(fit = readRDS(paste0("output/Fits/Sapling/",outloc,species,"/",species,"_fit.rds")),
                                 predictorstack = stacked,
                                 exclude = "yearmonth",
                                 BAspecies = species, # False if basal area is not included in vars, give species if ba is included
                                 get_dir = "data/Predictor_100m_Germany/",
                                 haconvert = TRUE,
                                 CI = TRUE, # for getting upper and lower confidence interval boundaries, Only use when link function is log
                                 cores = 20) 

  writeRaster(prediction, paste0("output/Predictions/Regeneration_",species,".tif"), overwrite=TRUE)
  }, FALSE)
}

print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))