# This script calibrates species-specific regeneration models while leaving out variable groups
# this variable group importance is tested only for final.species

start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", "DHARMa",
              "mgcViz", "tidyterra", "parallel", "cowplot",
              "modelr", "sf","blockCV", "automap", "stringr")
# sapply(packages, FUN = install.packages, character.only = T)
sapply(packages, FUN = library, character.only = T)


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
  
  
  # Leave category out model  ---------------------------------
for(leftout in unique(varimp$Category)){
  print(paste("left out:", leftout))
  
  # define fixed effects
  fixed = varimp %>% 
    filter(!(Category == leftout)) %>% # this category is ommited from fixed effects
    pull(varid)
  
  # create folder structure for saving
  dir.create(paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/"))
  for(species in species.final) {dir.create(paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species))}

  
  ## Start Model ---------------------------------------------------------------
  
  ### Iterate Species -----------------------------------------------------------
  for(species in species.final) {
    print(species)
    
    ### Built model variables ---------------------------------------------------
    mv <- modelVariables(species = species)
    mv %<>% # to avoid that more observations than the full set of variables
      select(all_of(c(resp.full, fixed.full, random.full, spatial.full, "plotid", "clusterid"))) %>% 
      drop_na()
    
    ### Run model ---------------------------------------------------------------
    try({
      fit <-  model.fit(resp = resp,
                        fixed = fixed,
                        random = random,
                        spatial = spatial,
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
                        CV = FALSE)
      saveRDS(fit, paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_fit.rds"))# save in directory
      
      ### Null model --------------------------------------------------------------
      print("Fit null model")
      fit.null <-  gam(formula = as.formula(paste(resp," ~ 1")), family = nb(), data = mv,
                       method = "REML",
                       select = FALSE)
      saveRDS(fit, paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_fitnull.rds"))# save in directory
      
      ### Calculate pseudo-R2 ------------------------------------------------------
      print("Calculate pseudo-R2")
      rsq.varimp = get_cohenrsq(fit = fit, 
                                fit_null = fit.null, 
                                test = mv,
                                train = mv,
                                resp = "count",
                                exclude = "yearmonth")
      saveRDS(rsq.varimp,  paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_rsq.rds"))
      print(paste("Done:", species))
    })
  }
  
  
  ## Model summary-------------------------------------------------------------
  # Save Model Summary in all species data frame
  df.allsp <- data.frame()
  
  for(species in species.final) {
    print(species)
    
    try({
      # set up output chart
      df.out <-  data.frame()
      
      # Load fit
      fit <- readRDS(paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_fit.rds"))
      
      # Get rsq
      rsq <- readRDS(paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/",species,"/",species,"_rsq.rds"))
      df.out[1, "rsq"] <- rsq$cohenrsq.train # test and train are the same
      
    })# close try
    
    # How many observations are used?
    df.out[1, "data.n"] <- length(fit$y)[1]
    
    # Save species
    df.out[1,"species"] <- species
    
    # Variable category
    df.out[1,"varcat"] <- leftout
    
    # save in all species df
    df.allsp <- bind_rows(df.allsp, df.out[1,])
  }
  
  # Save df
  write.csv2(df.allsp, paste0("output/Fits/Sapling/h50d7_Germany_varimp_",leftout,"/Sapling_model_summary.csv"), row.names = FALSE)
}


# Full model --------------------------------------------------------------
print("Full model")
start.full = Sys.time()

df.allsp.full <- data.frame()

for(species in species.final) {
  print(species)
  
  try({
    # set up output chart
    df.out <-  data.frame()
    
    # Data
    mv <- modelVariables(species = species)
    mv %<>% # to avoid that more observations than the full set of variables
      select(all_of(c(resp.full, fixed.full, random.full, spatial.full, "plotid", "clusterid"))) %>% 
      drop_na()
    
    # Load fit
    print("Load full model")
    fit <- readRDS(paste0("output/Fits/Sapling/h50d7_Germany/",species,"/",species,"_fit.rds"))
    
    
    fit.null <-  gam(formula = as.formula(paste(resp," ~ 1")), family = nb(), data = mv,
                     method = "REML",
                     select = FALSE)
    
    # Calculate pseudo R2
    print("Calculate pseudo-R2")
    rsq.full = get_cohenrsq(fit = fit, 
                            fit_null = fit.null, 
                            test = mv,
                            train = mv,
                            resp = "count",
                            exclude = "yearmonth")
    # Get rsq
    df.out[1, "rsq"] <- rsq.full$cohenrsq.train # test and train are the same
    
  })# close try
  
  # How many observations are used?
  df.out[1, "data.n"] <- length(fit$y)[1]
  
  # Species
  df.out[1,"species"] <- species
  
  # Variable category
  df.out[1,"varcat"] <- "full"
  
  # save in all species df
  df.allsp.full <- bind_rows(df.allsp.full, df.out[1,])
  print(paste("Done:",species, start.full-Sys.time()))
}

# Save df
write.csv2(df.allsp.full, paste0("output/Fits/Sapling/h50d7_Germany/Sapling_rsq_summary.csv"))

print(paste0("Full model done:",species, start.full-Sys.time()))

print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))