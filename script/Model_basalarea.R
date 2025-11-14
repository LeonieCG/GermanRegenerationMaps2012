# This script calibrates species-specific basal area models and checks model
# assumptions

start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", "DHARMa",
              "mgcViz", "tidyterra", "parallel", "cowplot",
              "modelr", "sf","blockCV", "automap")
# sapply(packages, FUN = install.packages, character.only = T)
sapply(packages, FUN = library, character.only = T)


# Load data ---------------------------------------------------------------

## Response Variables -------------------------------------------------------
WZP12 <-  readRDS("/bigdata/Inventories/DE BWI/Data/DE_BWI3_big_basalarea_wzp12.rds") 


## Explanatory variables ---------------------------------------------------
E <- readRDS("data/DE_BWI3_explvars_df.rds") # right one for the project
regeneration <- readRDS("data/DE_BWI3_regeneration_h50d7.rds") # is needed to built mv (has time in it...)

## Sapling species -----------------------------------------------------------------
source("script/Species_select.R")

for (i in species.vect) {dir.create(paste0("output/Fits/Basalarea/wzp12_interpol/",i))} # creates folder structure

## Source functions --------------------------------------------------------
source("script/Model_functions.R") # includes building Model variable table, Model function


# Model ---------------------------------------------------------------

## Iterate Species ---------------------------------------------------------------
for(species in species.vect) {
  print(species)

## Built model variables ---------------------------------------------------
  mv <- modelVariables(species = species)
  
  Data <- mv %>% 
    select(wzp12_ba_ha_species, x, y) %>% 
    drop_na()
  

## Run model ---------------------------------------------------------------
  try({
    start.model <- Sys.time()
    
    form_all <- as.formula("wzp12_ba_ha_species ~ s(x, y, bs = 'tp', k = 200)")
    print(form_all)
    
    fit <- bam(formula = form_all, family = tw(), data = Data
               , select = FALSE
               , method = "fREML" 
               , discrete=TRUE, nthreads=10) # speeds up calculation
    
    print(paste("Model took", Sys.time() - start.model, units(Sys.time()-start.model)))
    
    attr(fit, "species") <- attr(Data, "species")
    saveRDS(fit, paste0("output/Fits/Basalarea/wzp12_interpol/",species,"/Basalarea_",species,"_fit.rds"))
    })
}



# Model checks-------------------------------------------------------------

## Start Summary over all species ------------------------------------------------
# Save Model Summary in all species data frame
df.allsp <- data.frame()

# Save Model Summary in pdf
pdf(paste0("output/Fits/Basalarea/wzp12_interpol/Basalarea_DHARMaresidual.pdf"), width = 10, height = 10) # all graphs will be printed here
par(oma=c(1, 1, 1.5, 0.5),mfrow = c(2,1))

## Iterate Species ---------------------------------------------------------
for(species in species.vect) {
  print(species)
  try({
    fit  <-  readRDS(paste0("output/Fits/Basalarea/wzp12_interpol/",species,"/Basalarea_",species,"_fit.rds"))

    # set up output chart
    df.out <-  data.frame()

### Check DHARMa ------------------------------------------------------------
    sims <- simulateResiduals(fit, plot=T) #, exclude = paste0("s(", random, ")")) # for accounting for random effects add: exclude = paste0("s(", random, ")")
    
    mtext(paste0(species), cex = 1.5, outer = T, font = 3, side = 3, line = 0, adj = 1) # for plotting the species name

### Dipersion test ----------------------------------------------------------
    # for details see ?testDispersion
    disp <- testDispersion(sims) # tests under and overdispersion
    df.out[1,"disp.val"] <- round(disp$statistic, digits=3)
    df.out[1,"disp.p"] <- round(disp$p.value, digits=3)

### Zeroinflation test ------------------------------------------------------
    zeroinfl <- testZeroInflation(sims)
    df.out[1,"zeroinfl.val"] <- zeroinfl$statistic
    df.out[1,"zeroinfl.p"] <- zeroinfl$p.value
    
    mtext(paste0(species), cex = 1.5, outer = T, font = 3, side = 3, line = 0, adj = 1) # for plotting the species name
    
### Spatial autocorrelation test --------------------------------------------
    # due to the fact that all Ecken from the same Trakt (cluster)
    # have the same coordinates testSpatialAutocorrelation(sims) does not work.
    # Need for recalculation of residuals via clusterid.

    Data_group <- modelVariables(species = species) %>%
      select(all_of(c("clusterid", "wzp12_ba_ha_species","x","y"))) %>%
      drop_na() %>% # creates dataset with all important variables and deletes NAs
      group_by(clusterid) %>%
      arrange(x,y) %>%
      filter(row_number()==1) # keeps first row of every cluster group

    groupedSims <- recalculateResiduals(sims, group = Data_group$clusterid)# recalculated residuals

    spatautocorr <- testSpatialAutocorrelation(groupedSims, x = Data_group$x, y = Data_group$y)
    df.out[1,"spatautocorr.Iobserved"] <- spatautocorr$statistic[1]
    df.out[1,"spatautocorr.Iexpected"] <- spatautocorr$statistic[2]
    df.out[1,"spatautocorr.p"] <- spatautocorr$p.value
    
    mtext(paste0(species), cex = 1.5, outer = T, font = 3, side = 3, line = 0, adj = 1) # for plotting the species name

    })# close try

## Data stats --------------------------------------------------------------
    Data <- modelVariables(species = species) %>%
      select("wzp12_ba_ha_species","x","y") %>%
      drop_na()

    # How many plots are left when Nas dropped?
    df.out[1, "data.n"] <- dim(Data)[1]

    # How many plots are containing zeros?
    df.out[1, "data.resp.n0"] <- Data %>% select("wzp12_ba_ha_species") %>% summarise(sum(.==0))


## Save outputs -------------------------------------------------------------
    #df.out
    df.out[1,"species"] <- species
    write.csv(df.out, paste0("output/Fits/Basalarea/wzp12_interpol/",species,"/Basalarea_",species,"_fit.csv"), row.names = F)

    # save in all species df
    df.allsp <- bind_rows(df.allsp, df.out[1,])
}


## Save summary over all species -------------------------------------------
# Save df
write.csv2(df.allsp, "output/Fits/Basalarea/wzp12_interpol/Basalarea_model_summary.csv", row.names = FALSE)


#close pdf
dev.off()


print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))