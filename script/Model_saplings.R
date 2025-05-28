start.script <- Sys.time()

packages <- c("tidyr", "dplyr","magrittr", "terra", "mgcv", "gstat", "DHARMa",
              "mgcViz", "tidyterra", "parallel", "cowplot",
              "modelr", "sf","blockCV", "automap", "stringr")
# sapply(packages, FUN = install.packages, character.only = T)
sapply(packages, FUN = library, character.only = T)


# Choose output location --------------------------------------------------
outloc = "h50d7_Germany/" 

# Load data ---------------------------------------------------------------

## Response Variable -------------------------------------------------------
regeneration <- readRDS("data/DE_BWI3_regeneration_h50d7.rds")

## Explanatory variables ---------------------------------------------------
# Species independent
E <- readRDS("data/DE_BWI3_explvars_df.rds") 

# Species dependent
WZP12 <- readRDS("/bigdata/Inventories/DE BWI/Data/DE_BWI3_big_basalarea_wzp12.rds")

## Sapling species -----------------------------------------------------------------
source("script/Species_select.R")

## Source functions --------------------------------------------------------
source("script/Model_functions.R")  # includes building Model variable table, Model function, Cross validation

## Choose variables ----------------------------------------------------------
source("script/Model_vars.R")


# Start Model ---------------------------------------------------------------

## Iterate Species -----------------------------------------------------------
for(species in species.vect) {
  print(species)

## Built model variables ---------------------------------------------------
  mv <- modelVariables(species = species)

## Run model ---------------------------------------------------------------
  try({
  fit <-  model.fit(resp = resp,
                  fixed = fixed,
                  fixedfact = NULL,
                  random = random,
                  spatial = spatial,
                  offset = NULL,
                  exclude = "yearmonth",
                  fam = nb,
                  s.k = 10,
                  te.k= "c(25,50)",
                  ste.bs = "cs", # is faster than ts
                  select.var = FALSE,
                  bam = TRUE,
                  Data = mv,
                  CV = "blockcv",
                  blockcv.dr = 300000, # 300000 gives 11 blocks for whole germany
                  blockcv.k = 10)
    saveRDS(fit, paste0("output/Fits/Sapling/",outloc,species,"/",species,"_fit.rds"))
  })
}


# Model checks-------------------------------------------------------------

## Start Summary over all species ------------------------------------------------
# Save Model Summary in all species data frame
df.allsp <- data.frame()

# Save Model Summary in pdf
pdf(paste0("output/Fits/Sapling/",outloc,"/Sapling_DHARMaresidual.pdf"), width = 10, height = 10) # all graphs will be printed here
par(oma=c(1, 1, 1.5, 0.5),mfrow = c(2,1))

## Iterate Species ---------------------------------------------------------
for(species in species.vect) {
  print(species)

  try({
    # set up output chart
    df.out <-  data.frame()

    #Load fit
    fit <- readRDS(paste0("output/Fits/Sapling/",outloc,species,"/",species,"_fit.rds"))

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


## Check DHARMa ------------------------------------------------------------
    sims <- simulateResiduals(fit, plot=T) #, exclude = paste0("s(", random, ")")) # for accounting for random effects add: exclude = paste0("s(", random, ")")

    mtext(str_replace(species, "\\.", " "), cex = 1.5, outer = T, font = 3, side = 3, line = 0, adj = 1) # for plotting the species name

### Dispersion test ---------------------------------------------------------
    # for details see ?testDispersion
    disp <- testDispersion(sims) # tests under and overdispersion
    df.out[1,"disp.val"] <- round(disp$statistic, digits=3)
    df.out[1,"disp.p"] <- round(disp$p.value, digits=3)

### Zeroinflation test-------------------------------------------------------------------------
    zeroinfl <- testZeroInflation(sims)
    df.out[1,"zeroinfl.val"] <- zeroinfl$statistic
    df.out[1,"zeroinfl.p"] <- zeroinfl$p.value

    mtext(str_replace(species, "\\.", " "), cex = 1.5, outer = T, font = 3, side = 3, line = 0, adj = 1) # for plotting the species name

### Spatial autocorrelation test -------------------------------------------------------------------------
    # due to the fact that all Ecken (plots) from the same Trakt (cluster)
    # have the same coordinates testSpatialAutocorrelation(sims) does not work.
    # Need for recalculation of residuals via clusterid.

    Data_group <- modelVariables(species = species) %>%
      select(all_of(c("clusterid", resp, fixed, random, spatial))) %>%
      drop_na() %>% # creates dataset with all important variables and deletes NAs
      group_by(clusterid) %>%
      arrange(x,y) %>%
      filter(row_number()==1) # keeps first row of every cluster group

    groupedSims <- recalculateResiduals(sims, group = Data_group$clusterid)# recalculated residuals

    spatautocorr <- testSpatialAutocorrelation(groupedSims, x = Data_group$x, y = Data_group$y)# plot = FALSE
    df.out[1,"spatautocorr.Iobserved"] <- spatautocorr$statistic[1]
    df.out[1,"spatautocorr.Iexpected"] <- spatautocorr$statistic[2]
    df.out[1,"spatautocorr.p"] <- spatautocorr$p.value

    mtext(str_replace(species, "\\.", " "), cex = 1.5, outer = T, font = 3, side = 3, line = 0, adj = 1) # for plotting the species name
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
  write.csv(df.out, paste0("output/Fits/Sapling/",outloc,species,"/",species,"_fit.csv"), row.names = F)

  # save in all species df
  df.allsp <- bind_rows(df.allsp, df.out[1,])

}


## Save summary over all species -------------------------------------------
# Save df
write.csv2(df.allsp, paste0("output/Fits/Sapling/",outloc,"Sapling_model_summary.csv"), row.names = FALSE)

#close pdf
dev.off()


print(paste("Script took",Sys.time()-start.script, units(Sys.time()-start.script)))