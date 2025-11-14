# This script predicts species-specific basal area for Germany

start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", "DHARMa",
              "mgcViz", "tidyterra", "parallel")
sapply(packages, FUN = library, character.only = T)

# Load data ---------------------------------------------------------------
# Sapling species
source("script/Species_select.R")


# Source functions --------------------------------------------------------
source("script/Prediction_function.R")


# Choose variables ----------------------------------------------------------
source("script/Model_vars.R")


# Predictor stack ---------------------------------------------------------
stacked <- c(rast("data/Predictor_100m_Germany/x.tif"),
             rast("data/Predictor_100m_Germany/y.tif"))

# Run prediction ----------------------------------------------------------
for(species in species.vect) {
  print(species)
  try({
    prediction = predict.fit(fit = readRDS(paste0("output/Fits/Basalarea/wzp12_interpol/",species,"/Basalarea_",species, "_fit.rds")),
                             predictorstack = stacked,
                             exclude = NULL ,#"yearmonth",
                             BAspecies = NULL,
                             get_dir = NULL,
                             haconvert = FALSE, # converts to hectare
                             CI = FALSE, # Only use when link function is log
                             cores = 20) # for parallel
    writeRaster(prediction, paste0("data/Predictor_100m_Germany/wzp12_ba_ha_species_interpol_",species,".tif"), overwrite=TRUE)
  }, FALSE)
 }


print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))
