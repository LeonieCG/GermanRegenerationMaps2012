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
    prediction = predict.fit(fit = readRDS(paste0("output/Fits/Basalarea/wzp12/",species,"/Basalarea_",species, "_fit.rds")),
                             predictorstack = stacked,
                             exclude = NULL ,
                             BAspecies = NULL,
                             get_dir = NULL,
                             haconvert = FALSE, # converts to hectare
                             CI = FALSE, # Only use when link function is log
                             cores = 20) # for parallel
    writeRaster(prediction, paste0("data/Predictor_100m_Germany/wzp12_ba_ha_species_",species,".tif"), overwrite=TRUE)
  }, FALSE)
 }



# Retrieve predicted basal area values at plot location ------------------------------------
loc <- read.csv("data/bwi_wzp12_ba_ha_cells.csv") %>% # cell ids extracted by Thünen
  mutate(plotid = paste0("DE_BWI_",tnr,".",enr)) %>% 
  select(plotid, cell)

for(species in species.vect){
  ba  <- rast(paste0("data/Predictor_100m_Germany/wzp12_ba_ha_species_",species,".tif"))
  vals <- extract(ba, loc$cell)# extract basal area
  loc <- cbind(loc, vals)# combine
}

wzp12_interpol <- loc %>% 
  select(-cell) %>% 
  pivot_longer(cols = 2:dim(.)[2], names_to = "tax", values_to = "wzp12_ba_ha_species") %>% 
  group_by(tax) %>% 
  nest()

saveRDS(wzp12_interpol, "data/DE_BWI3_big_basalarea_wzp12_interpol.rds")

print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))
