# This script predicts species-specific basal area for Germany

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

# due to licencing, not all environmental predictors could be published with this study 
# the following is a selection, model outcomes can thereby vary highly from the study!
resp = "count"
fixed = c("tPeriodic2020_forcli", "tSeas2020_forcli","tMinColdMonth2020_forcli","tRangeDay2020_forcli","tRangeAn2020_forcli" #T microclimate
          , "tPeriodic2010_chelsa", "tSeas2010_chelsa","tMinColdMonth2010_chelsa","tRangeDay2010_chelsa","tRangeAn2010_chelsa" #T macroclimate
          , "precPeriodic2010_chelsa", "precSeas2010_chelsa" # prec
          , "wwpi_cop" # Water prob index, Anoxy indicator, water bodies and flooded plains
          , "tcd_cop" #tree cover density
          , "alt", "northexp", "eastexp" # terrain vars
          , "wzp12_ba_ha_species" #Basal area of respective old trees
)
random = c("yearmonth" # to account for changes within sampling period
           , "blname" # Bundesland
)
spatial = c("x","y")

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
