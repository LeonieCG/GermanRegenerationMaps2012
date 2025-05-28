start <- Sys.time()

# Packages ----------------------------------------------------------------
packages <- c("tidyr", "dplyr", "stringr","magrittr","terra","sf", "tidyterra", "Predictors", "Hmisc",
              "ggplot2")
# sapply(packages, FUN = install.packages)
sapply(packages, FUN = library, character.only = T)


# Data --------------------------------------------------------------------
## Species -----------------------------------------------------------------
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")

species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species


# All species that have good models and a cultivation risk map
species.final.cult <- 
  species.tab %>% 
  filter(name.id %in% species.final) %>% 
  filter(!is.na(name.cultrisk)) %>% 
  dplyr::pull(name.id)


## Forest area - Bavaria ---------------------------------------------------
# in order to have the same coordinate system rather project by.f and regeneration maps, since reprojecting leads to new res and to recalculation of values,
# and the Anbaurisikokarten do have categorical values not continuous! Project to get same extent!
bay.f <- rast("data/ForestArea_bavaria_100m.tif") %>% 
  rename(bay_f = names(.))

bay.f %<>% 
  project(., rast("data/Anbaurisikokarten_basis01_lwf/basisLite_abr_final_Tcord2100.tif"))


# Prepare data -------------------------------------------------------------------------
df_total = data.frame()

for (species in species.final.cult){
print(species)

## Cultivation risk -Bavaria -----------------------------------------------
# Cultivation risk categories
# 11-13 = sehr geringes Risiko = very low risk
# 21-23 = geringes Risiko = low risk
# 31-33 = erhöhtes Risiko = increased risk
# 41-43 = hohes Risiko = high risk
# 51-53 = sehr hohes Risiko = very high risk

if(species == "Tilia"){ # Tilia has two species combined in BWI, combine cultrisk maps
  tcord <- rast("data/Anbaurisikokarten_basis01_lwf/basisLite_abr_final_Tcord2100.tif") %>%
    rename(cultrisk_val = names(.)) %>%
    mutate(cultrisk_val = as.numeric(cut(cultrisk_val, breaks = c(11,13,23,33,43,53), labels =  c("1", "2", "3", "4", "5"), include.lowest = TRUE)))

  tplat  <- rast("data/Anbaurisikokarten_basis01_lwf/basisLite_abr_final_Tplat2100.tif") %>%
    rename(cultrisk_val = names(.)) %>%
    mutate(cultrisk_val = as.numeric(cut(cultrisk_val, breaks = c(11,13,23,33,43,53), labels =  c("1", "2", "3", "4", "5"), include.lowest = TRUE)))

  risk <- round(mean(tcord,tplat)) %>% # conservative mean of both
    mutate(cultrisk_en = cut(cultrisk_val, breaks = c(1,3,5), labels =  c("lower", "higher"), include.lowest = TRUE))%>%
    select(cultrisk_en)

  } else {
  risk <- rast(paste0("data/Anbaurisikokarten_basis01_lwf/basisLite_abr_final_",species.tab[species.tab$name.id == species, "name.cultrisk"],"2100.tif"))

  risk %<>%
    rename(cultrisk_val = names(.)) %>%
    # mutate(cultrisk_en = cut(cultrisk_val, breaks = c(11,13,23,33,43,53), labels =  c("very_low", "low", "increased", "high", "very_high"), include.lowest = TRUE))
    mutate(cultrisk_en = cut(cultrisk_val, breaks = c(11,33,53), labels =  c("lower", "higher"), include.lowest = TRUE)) %>%
    select(cultrisk_en)
}

## Regeneration Prediction -------------------------------------------------
reg <- rast(paste0("output/Predictions/Regeneration_",species,".tif")) %>%
  rename(count_ha = paste0(species)) %>%
  select(count_ha)

# Coordinate system
# rather reproject regeneration map since it is continuous and by recalculation resolution a mean value can be built,
# not possible for factor e.g. cultivation risk
reg %<>%
  project(., risk)

## Merge all data and create data.frame ------------------
# Merge Forest area, cultivation risk and regeneration
regrisk <- c(bay.f,risk,reg)

regrisk.df <- regrisk %>%
  filter(bay_f == 1) %>% # select just forest cells
  as.data.frame(.,cells=TRUE) %>% # turn into data frame
  mutate(species = species) # change species cell name to "species"

# add species df to df with including all species
df <- regrisk.df
df_total <- rbind(df_total,df)
}

# save preliminary results
saveRDS(df_total, "output/Cultivationrisk/df_total_cache.rds")
df_total <- readRDS("output/Cultivationrisk/df_total_cache.rds")


# Complete for cell, cultrisk and species calculate ------------------------------------------------------------------
# TEST data
# df_total <- data.frame(cell =c(1,2,3,1,2,3,1,2,3),
#                        species=c("FS","FS","FS","AA","AA","AA","AG","AG","AG"),
#                        cultrisk_en=as.factor(c("lower","lower","lower","higher","higher","higher","NA","higher","lower")),
#                        count_ha=c(1,3,2,4,20,NA,23,NA,NA))

df_total_sum <-
  df_total %>% 
  group_by(cell, cultrisk_en) %>% 
  summarise(count_ha = sum(count_ha, na.rm= F), .groups = "drop") %>% # count_ha per cell and cultrisk
  group_by(cell) %>% 
  mutate(count_ha_cell = sum(count_ha,na.rm = F),
         count_percent = count_ha/count_ha_cell*100) %>% #percentage saplings per cult risk category and cell
  complete(cultrisk_en,
           fill=list(count_ha = NA, count_ha_cell = NA, count_percent = NA)) %>% # complete for each cell and cultrisk category
  ungroup()

# save preliminary results
saveRDS(df_total_sum, "output/Cultivationrisk/df_total_sum_cache.rds")
df_total_sum <- readRDS("output/Cultivationrisk/df_total_sum_cache.rds")


# Sapling cultivation risk spatRaster -------------------------------------
# change from long to wide table, for each cultrisk category own column with count_percent as cell values
scr <- df_total_sum %>% 
  select(c(cell, cultrisk_en, count_percent)) %>% 
  pivot_wider(., names_from = cultrisk_en, values_from = count_percent)

# get coordinates and cell number from forest area of bavaria
rast_df <- bay.f %>%
  as.data.frame(.,cells=TRUE, xy=T) %>% 
  select(cell,x,y)

# merge forest area coodrinates with cultrisk via cell number
scr <- merge(rast_df, scr, by = "cell", all = T) %>%
  select(-c(cell))

# turn data.frame into a spatRaster via coordinates
scr_rast <- as_spatraster(scr, xycols = 1:2, crs = crs(bay.f))
plot(scr_rast)

writeRaster(scr_rast, "output/Cultivationrisk/Regeneration_cultivationrisk.tif", overwrite = TRUE)


end <- Sys.time()
print(end-start)
