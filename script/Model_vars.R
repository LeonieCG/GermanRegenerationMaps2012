## Specify explanatory and prediction variables
resp = "count"
fixed = c("tPeriodic2020_forcli", "tSeas2020_forcli","tMinColdMonth2020_forcli","tRangeDay2020_forcli","tRangeAn2020_forcli" #T microclimate
          , "tPeriodic2010_chelsa", "tSeas2010_chelsa","tMinColdMonth2010_chelsa","tRangeDay2010_chelsa","tRangeAn2010_chelsa" #T macroclimate
          , "precPeriodic2010_chelsa", "precSeas2010_chelsa" # prec
          , "wwpi_cop" # Water prob index, Anoxy indicator, water bodies and flooded plains
          , "ai_cgiar"# aridity index
          , "clay_esdact", "silt_esdact", "sand_esdact","bulkd_esdact", "coarse_esdact" #Soil properties
          , "awc_esdact" # "Soil water balance"
          , "nfkwe_bgr" # nFK im effektiven Wurzelraum
          , "octop_esdacoc" # OC, is there peat?
          , "k_esdacc", "n_esdacc", "p_esdacc" # Soil Nutrients
          , "phCaCl_esdacc","cn_esdacc","cec_esdacc", "caco3_esdacc" #Soil chemistry and nutrition availability
          , "nTotalFlux_pineti", "nh4TotalFlux_pineti", "no3TotalFlux_pineti" # nitrogen, nh4, no3 immission
          , "pet_cgiar" # potential evapotranspiration
          , "gdd0_aclim", "cwbGdd0_aclim", "cwbYear_aclim" #cwb=climatic water balance, cwbYear_aclim seems more differentiated
          , "tcd_cop" #tree cover density
          , "alt", "northexp", "eastexp" # terrain vars
          , "wzp12_ba_ha_species" #Basal area of respective old trees
          )
random = c("yearmonth" # to account for changes within sampling period
           , "blname" # Bundesland
           )

spatial = c("x","y")
