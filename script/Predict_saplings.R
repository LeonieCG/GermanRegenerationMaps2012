# This script predicts species-specific regeneration for Germany

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