start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", "DHARMa",
              "mgcViz", "tidyterra", "parallel")
sapply(packages, FUN = library, character.only = T)

# Load data ---------------------------------------------------------------
# Sapling species
species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species

# Source functions --------------------------------------------------------
source("script/Prediction_function.R")


# Choose variables ----------------------------------------------------------
source("script/Model_vars.R")
resp = "wzp12_ba_ha_species" 
fixed = fixed[! fixed %in% c("wzp12_ba_ha_species")]# remove wzp12 because it is not needed


# Predictor stack ---------------------------------------------------------
stacked <-  predict.stack(vars =  c(fixed, spatial, random),
                          exclude = "yearmonth",
                          get_dir = "data/Predictor_100m_Germany/") 



# Run prediction ----------------------------------------------------------
for(species in species.final) {
  print(species)
  try({
    prediction = predict.fit(fit = readRDS(paste0("output/Fits/Basalarea/wzp12/",species,"/Basalarea_",species, "_fit.rds")),
                             predictorstack = stacked,
                             exclude = "yearmonth",
                             BAspecies = NULL,
                             get_dir = NULL,
                             haconvert = FALSE, # converts to hectare
                             CI = FALSE, # Only use when link function is log
                             cores = 20) # for parallel
    writeRaster(prediction, paste0("data/Predictor_100m_Germany/wzp12_ba_ha_species_",species,".tif"), overwrite=TRUE)
  }, FALSE)
 }


print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))
