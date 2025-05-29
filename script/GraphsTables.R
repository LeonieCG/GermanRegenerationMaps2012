# Graphics, tables and summary statistics for publication

# Packages ----------------------------------------------------------------
packages <- c("dplyr","tidyr","magrittr", "terra", "tidyterra","viridis", "ggplot2", 
              "stringr", "cowplot", "patchwork", "ggExtra")
sapply(packages, FUN = library, character.only = T)


# data --------------------------------------------------------------------
dt.f <- rast("data/ForestArea_germany_100m.tif")
bay.f <- rast("data/ForestArea_bavaria_100m.tif")

germany <- vect("data/DE_Verwaltungsgebiete5000/vg5000_ebenen_1231/VG5000_KRS.shp") %>% 
  aggregate() %>% 
  project(., dt.f)
bavaria <- vect("data/DE_Verwaltungsgebiete5000/vg5000_ebenen_1231/VG5000_LAN.shp") %>% 
  filter(GEN == "Bayern") %>% 
  project(.,dt.f)

# colors
sunset = colorRampPalette(c("#FFEC9DFF", "#F2AF4AFF", "#EB7F54FF", "#C36377FF", "#61599DFF", "#1D457F", "#191F40FF", "black"))
cult.col = colorRampPalette(c("grey90", "#FCFD8F","#F3CE65","#EB9F3C","#9A3F07"))
div = colorRampPalette(c("grey90","#FFEC9DFF", "#F2AF4AFF", "#EB7F54FF", "#9A3F07"))


# REGENERATION DISTRIBUTION ------------------------------------------
# Visualising predicted species distribution Maps

## Distribution Species ----------------------------------------------------------------
sunset = colorRampPalette(c("#FFEC9DFF", "#F2AF4AFF", "#EB7F54FF", "#C36377FF", "#61599DFF", "#1D457F", "#191F40FF", "black"))

### Plot --------------------------------------------------------------------
sapling.map <- function(species.vect=species.vect, source, scale, scale.plot, max.count){
  for(species in species.vect){
    print(species)
    # load
    pred <- rast(paste0("output/Predictions/Regeneration_",species,".tif")) %>% 
      rename(count_ha = paste0(species)) %>% 
      select(count_ha)
      # Sapling
    if (scale == "Bavaria") {
      pred %<>% terra::crop(., scale.plot, mask = T)}
    # plot
    p <-
      ggplot()+
      theme_void()+
      geom_spatraster(data = pred) +
      scale_fill_gradientn(
        "Regeneration\ndensity [ha\u207B\u00B9]",
        na.value = "transparent",
        colors = sunset(7),
        space = "Lab",
        trans = "log1p",
        breaks = c(0,10,100,1000,10000, 100000, 1e6,1e7),#Sapling
        labels = scales::comma(c(0,10,100,1000,10000, 100000,1e6,1e7)),#Sapling
        limits = c(0, max.count))+
      guides(fill = guide_colourbar(barwidth = 1))+
      theme(legend.title = element_text(size = 8),
            legend.text = element_text(size = 7))+
      geom_spatvector(data = scale.plot, fill = "transparent", colour = "black", linewidth = 0.1)+
      annotate("text", x = -Inf, y = Inf, hjust=-0.1, vjust = 1, size = 3, label = paste(species.tab[species.tab$name.id==species,]$name.scient), fontface = 'bold.italic')
    
    assign(paste0("p.",species), p , envir = .GlobalEnv)
      
    # save
    ggsave(plot = p,
           filename = paste0("output/Graphs/Regeneration_",scale,"_",species,".png"),
           height = 8, width = 8, units = "cm", dpi = 900,
           bg = "white",
           device=grDevices::png)
  }
}

species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")

# choosing max.count
regstack <- rast(paste0("output/Predictions/Regeneration_",species.final[1],".tif")) %>%
  select(paste0(species.final[1]))

for (i in 2:length(species.final)){
  stackadd <- rast(paste0("output/Predictions/Regeneration_",species.final[i],".tif")) %>%
    select(paste0(species.final[i]))
  regstack <- c(regstack, stackadd)
  rm(stackadd)
}

global(regstack, "max",na.rm=T)

# 99% Quantile
quantile(terra::values(regstack), 0.99, na.rm = TRUE)# of all regeneration maps

regstack.rest <- regstack %>%
  select(species.final[!species.final %in% c("Picea.abies","Abies.alba", "Fagus.sylvatica")])
quantile(terra::values(regstack.rest), 0.99, na.rm = TRUE) # what cut off for the color scale makes sense earth engine question


## Final species -----------------------------------------------------------
#sapling.map(scale = "Bavaria", scale.plot = bavaria, addcoord = T, max.count = 21000)
sapling.map(species.vect = species.final, 
            scale = "Germany", scale.plot = germany, max.count = 1186) #cut off at 99% quantile

paste0("p.",sort(species.final[!species.final %in% c("Picea.abies","Abies.alba", "Fagus.sylvatica")]), collapse = " + ")

p.Acer.campestre + p.Acer.platanoides + p.Alnus.glutinosa + p.Alnus.incana + p.Betula.pubescens +
  p.Carpinus.betulus + p.Castanea.sativa + p.Fraxinus.excelsior + p.Larix.kaempferi +
  p.Pinus.sylvestris + p.Populus.nigra + p.Prunus.avium + p.Prunus.serotina + p.Pseudotsuga.menziesii +
  p.Quercus.robur + p.Quercus.rubra + p.Robinia.pseudoacacia + p.Sorbus.aria + p.Tilia +
  plot_layout(ncol = 4, nrow = 5, guides = 'collect')

ggsave(filename = "output/Graphs/Regeneration_Germany_rest.png",
       height = 28, width = 20, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)

## Selected species --------------------------------------------------------
sapling.map(species.vect = c("Picea.abies","Abies.alba", "Fagus.sylvatica"), # max count is 44073.406 of FS
            scale = "Germany", scale.plot = germany, max.count = 45000) 

p.Fagus.sylvatica + p.Picea.abies + p.Abies.alba + plot_layout(guides = 'collect')

ggsave(filename = "output/Graphs/Regeneration_Germany_AA_FS_PA.png",
       height = 7, width = 17, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)


# TOTAL DENSITY --------------------------------------------------------------------
# Calculation -------------------------------------------------------------
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")
species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species
# species.final <- "Abies.alba" # for test reasons
# species = "Abies.alba"

regstack <- rast(paste0("output/Predictions/Regeneration_",species.final[1],".tif")) %>%
  select(paste0(species.final[1]))

for (i in 2:length(species.final)){
  stackadd <- rast(paste0("output/Predictions/Regeneration_",species.final[i],".tif")) %>%
    select(paste0(species.final[i]))
  regstack <- c(regstack, stackadd)
  rm(stackadd)
}

regtot = sum(regstack)
names(regtot) <- "count_tot_ha"

regtot$class <- regtot %>% classify(c(0,1000,2000,Inf),right = F)

regtot.df = as.data.frame(regtot)

saveRDS(regtot.df, "output/Graphs/Density.rds")
writeRaster(regtot, "output/Graphs/Density.tif", overwrite=TRUE)


regtot.df <- readRDS("output/Graphs/Density.rds") 
regtot <- rast("output/Graphs/Density.tif")

## Plot --------------------------------------------------------------------

### Class -------------------------------------------------------------------
p.tot.class <- 
  ggplot() + 
  theme_void() +
  geom_spatraster(data = regtot$class)+
  scale_fill_manual(
    "Regeneration\ndensity [ha\u207B\u00B9]",
    na.translate = FALSE,
    labels = c("0-1,000", "1,000-2,000", "\u22652,000"),
    values = c("#9A3F07","#F2AF4AFF","grey90"),
    guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = germany, fill = "transparent", colour = "black", linewidth = 0.1)

# Save
ggsave(plot = p.tot.class,
       filename = paste0("output/Graphs/Density.png"),
       height = 8, width = 8, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)


### Continuous --------------------------------------------------------------
p.tot.cont <- 
  ggplot()+
  theme_void()+
  geom_spatraster(data = regtot$count_tot_ha) +
  scale_fill_gradientn(
    "Regeneration\ndensity [ha\u207B\u00B9]",
    na.value = "transparent",
    colors = sunset(7),
    space = "Lab",
    trans = "log1p",
    breaks = c(0,10,100,1000,10000, 100000, 1e6,1e7),#Sapling
    labels = scales::trans_format("log10", scales::math_format(10^.x)),#Sapling
    limits = c(0, global(regtot$count_tot_ha, max, na.rm=T)$max))+
  guides(fill = guide_colourbar(barwidth = 1))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  geom_spatvector(data = germany, fill = "transparent", colour = "black", linewidth = 0.1)

p.tot.cont.hist <-
  ggplot(regtot.df, aes(x = count_tot_ha, fill = ..x..)) +
  theme_minimal()+
  geom_histogram(bins = 40,
                 boundary = 0,
                 color="black",
                 linewidth = 0.1)+
  xlab("Regeneration density [ha\u207B\u00B9]") +
  ylab("Area [ha]") +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_x_continuous(trans = "log1p",
                     limits = c(0, NA),
                     breaks = c(0,10,100,1000,10000, 100000, 1e6,1e7),#Sapling
                     labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_fill_gradientn(colors = sunset(7),
                       trans = "log1p",
                       limits = c(0, global(regtot$count_tot_ha, max, na.rm=T)$max),
                       breaks = c(0,10,100,1000,10000, 100000, 1e6,1e7),#Sapling
                       labels = scales::trans_format("log10", scales::math_format(10^.x)),
                       guide = "none")
layout <- c(
  area(t = 0, b = 10, l = 0,  r = 12),
  area(t = 0, b = 10, l = 13, r = 20))

free(p.tot.cont) + p.tot.cont.hist + plot_layout(guides = 'collect') +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "output/Graphs/S_Density.png",
       height = 8, width = 17, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)


## Total count stats -------------------------------------------------------
# Mean count per hectare
mean(regtot.df$count_tot_ha,na.rm=T)

# how much area is affected by < 1000 , 1000-2000, 20000
length(regtot.df[regtot.df$count_tot_ha < 1000,])/length(regtot.df$count_tot_ha)*100
length(regtot.df[regtot.df$count_tot_ha >= 1000 & regtot.df$count_tot_ha < 2000,])/length(regtot.df$count_tot_ha)*100
length(regtot.df[regtot.df$count_tot_ha >= 2000,])/length(regtot.df$count_tot_ha)*100

## Table mean density per species --------------------------------------------------
reg.stats = global(c(regstack,regtot), c("mean","sd"), na.rm=TRUE)
write.csv2(reg.stats, "output/Graphs/Regeneration_density_stats.csv")


# SPECIES RICHNESS ---------------------------------------------------------------
# species richness rule BaySF Waldbauhandbuch Baumartenwahl 2020
# 4 species per hectare for species with minimum abundances of >=5 %
# 3 species Lindner et al 2020 (policy brief)

regtot = sum(regstack)
names(regtot) <- "count_tot_ha"
regstack.regtot.rast = c(regstack, regtot)
regstack.regtot.df = as.data.frame(regstack.regtot.rast, xy=T)

sprich.df <- regstack.regtot.df%>%
  mutate(across(.cols = -c(x,y,count_tot_ha), .fns = function(x) x/count_tot_ha)) %>%
  select(-c(count_tot_ha)) %>%
  mutate(across(.cols = -c(x,y), .fns = function(x) ifelse(x < 0.05, 0, 1))) %>%
  mutate(sprich = rowSums(across(.cols = -c(x,y)))) %>%
  select(c(x,y,sprich))

sprich.rast <-  as_spatraster(sprich.df, xycols = 1:2, crs = crs(regstack))

sprich.rast$class<- sprich.rast %>%
  classify(c(0,2,4,Inf)) # 0-2 low, 3-4 intermediate, >=5 enough

saveRDS(as.data.frame(sprich.rast), "output/Graphs/Speciesrichness.rds")
writeRaster(sprich.rast, "output/Graphs/Speciesrichness.tif", overwrite=TRUE)

sprich.df <- readRDS("output/Graphs/Speciesrichness.rds")
sprich.rast <- rast("output/Graphs/Speciesrichness.tif")
  
## Plot --------------------------------------------------------------------

### Class -------------------------------------------------------------------

p.div <- 
  ggplot() +
  theme_void() +
  geom_spatraster(data = sprich.rast$class) +
  scale_fill_manual(
    "",
    na.value = "transparent",
    # labels = c("0-3", "4", "5-11",""),
    values = c("#9A3F07","#F2AF4AFF","grey90"),
    guide = 'none') +
  geom_spatvector(data = germany, fill = "transparent", colour = "black", linewidth = 0.1)

p.div.hist <-
  ggplot(sprich.df, aes(x = sprich, fill =..x..)) +
  theme_classic() +
  xlab("Species richness") +
  ylab("Area [10\u2076 ha]") +
  geom_histogram(binwidth = 1 ,
                 boundary=-0.5,
                 color="black") +
  scale_x_continuous(breaks = 1:global(sprich.rast$sprich, max, na.rm=T)$max) +
  scale_fill_gradientn("",
    na.value = "transparent",
    colors = c("#9A3F07","#9A3F07","#F2AF4AFF","#F2AF4AFF","grey90","grey90"), 
    values = scales::rescale(c(0,1.99,2,3.99,4,11)),
    guide = 'none') +
  scale_y_continuous(limits=c(0,35e5),
                     breaks = c(0,1e6,2e6,3e6),
                     labels = scales::comma(c(0,1,2,3))) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  annotate("segment", x = 2.5, xend = 0.5, y = 35e5, yend = 35e5,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F") +
  annotate("text", x = 2.7, y = 35e5, 
           hjust = 0, vjust = 0.5,
           label = paste(strwrap(paste0(round(table(sprich.df$sprich <= 2)[2]/dim(sprich.df)[1]*100, digits = 1),"% of the forest area has a tree species richness \u22642"), 30), collapse = "\n"),
           colour = "#1D457F",
           size = 2)

layout <- c(
      area(t = 0, b = 10, l = 0,  r = 12),
      area(t = 0, b = 10, l = 13, r = 20))

free(p.div) + p.div.hist + 
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")
  
ggsave(filename = "output/Graphs/Speciesrichness.png",
       height = 8, width = 11, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)

### Continuous --------------------------------------------------------------
p.div.cont <- 
  ggplot()+
  theme_void()+
  geom_spatraster(data = sprich.rast$sprich) +
  scale_fill_gradientn(
    "Species richness",
    na.value = "transparent",
    colors = sunset(7),
    space = "Lab",
     breaks = c(0,2,4,6,8,10),#Sapling
    # labels = scales::trans_format("log10", scales::math_format(10^.x)),#Sapling
    limits = c(0, global(sprich.rast$sprich, max, na.rm=T)$max))+
  guides(fill = guide_colourbar(barwidth = 1))+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))+
  geom_spatvector(data = germany, fill = "transparent", colour = "black", linewidth = 0.1)

p.div.cont.hist <-
  ggplot(sprich.df, aes(x = sprich, fill = ..x..)) +
  theme_minimal()+
  geom_histogram(binwidth = 1,
                 color="black",
                 linewidth = 0.1)+
  xlab("Species richness") +
  ylab("Area [ha]") +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_x_continuous(breaks = 0:12,
                     labels = 0:12) +
  scale_fill_gradientn(colors = sunset(7),
                       limits = c(0,NA),
                       breaks= 0:11,
                       guide = "none")

layout <- c(
  area(t = 0, b = 10, l = 0,  r = 12),
  area(t = 0, b = 10, l = 13, r = 20))

free(p.div.cont) + p.div.cont.hist + plot_layout(guides = 'collect') +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

ggsave(filename = "output/Graphs/S_Speciesrichness.png",
       height = 8, width = 17, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)



## Stats -------------------------------------------------------------------
# Average
mean(sprich.df$sprich)

# how much of the forest has a sp richness <=2?
table(sprich.df$sprich <=2)[2]/dim(sprich.df)[1]*100

# how much of the forest has a sp richness >= 5?
table(sprich.df$sprich >= 5)[2]/dim(sprich.df)[1]*100

# how much of the forest has a sp richness 3-4?
(table(sprich.df$sprich == 3)[2]+ table(sprich.df$sprich == 4)[2])/dim(sprich.df)[1]*100

# how much of the forest has a sp richness > 4?
table(sprich.df$sprich >= 3)[2]/dim(sprich.df)[1]*100


# CULTIVATION RISK ---------------------------------------------------------
## Species -----------------------------------------------------------------
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")

species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species
# species.final <- "Abies.alba" # for test reasons
# species = "Abies.alba"

cult.col = colorRampPalette(c("grey90", "#FCFD8F","#F3CE65","#EB9F3C","#9A3F07"))

sap.risk <- rast("output/Cultivationrisk/Regeneration_cultivationrisk.tif") %>% 
  select(higher)
crs(sap.risk)==crs(bavaria)


## Plot ----------------------------------------------------------------
### Map ----------------------------------------------------------------
p.cult <-  
  ggplot()+
  theme_void()+
  geom_spatvector(data = bavaria, fill = "white", colour = "transparent")+
  geom_spatraster(data = sap.risk) +
  scale_fill_gradientn("",
                       colours = cult.col(50),
                       na.value = "transparent",
                       limits = c(0, 100),
                       guide = "none") +
  geom_spatvector(data = bavaria, fill = "transparent", colour = "black", linewidth = 0.1)


### Map histogram ------------------------------------------------------
sap.risk.df.nona <- readRDS("output/Cultivationrisk/df_total_sum_cache.rds") %>% 
  select(c(cell, cultrisk_en, count_percent)) %>% 
  pivot_wider(., names_from = cultrisk_en, values_from = count_percent) %>% 
  select(c(higher,cell)) %>% 
  drop_na()

# Frequency
h.cult <-
  ggplot(sap.risk.df.nona, aes(x = higher)) +
  theme_classic() +
  geom_histogram(binwidth = 5,
                 boundary = 0,#-0.5
                 fill=cult.col(20), color="black") +
  xlab("Proportion of regeneration at\nhigh cultivation risk [%]") +
  ylab("Area [10\u00B3 ha]")+
  scale_y_continuous(limits=c(0,69e4),
                     breaks = c(0,2e5,4e5,6e5),
                     labels = c(0,200,400,600)) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  annotate("segment", x = round(median(sap.risk.df.nona$higher, na.rm =T), digits=0), xend = 0, y = 65e4, yend = 65e4,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F") +
  annotate("text", x = round(median(sap.risk.df.nona$higher, na.rm =T))+ 2, y = 65e4, 
           hjust = 0, vjust = 0.5,
           label = paste(strwrap(paste0("50% of forest area has a cultivation risk <", round(median(sap.risk.df.nona$higher, na.rm =T), digits=1),"%"), 20), collapse = "\n"),
           colour = "#1D457F",
           size = 2) +
  annotate("segment", x = 100, xend = 75, y = 2e5, yend = 2e5,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F")+
  annotate("text", x = 100, y = 2.5e5, 
           hjust = 1, vjust = -0.1,
           label = paste(strwrap(paste0(round(table(sap.risk.df.nona$higher >= 75)[2]/dim(sap.risk.df.nona)[1]*100, digits=1),"% forest area has a high proportion of regeneration at risk (\u226575%)"), 25), collapse = "\n"),
           colour = "#1D457F",
           size = 2)

table(sap.risk.df.nona$higher >= 75)[2]

layout <- c(
  area(t = 0, b = 10, l = 0,  r = 13),
  area(t = 0, b = 10, l = 14, r = 20))

free(p.cult) + h.cult +
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "A")

ggsave(filename = paste0("output/Graphs/Cultivationrisk_regeneration.png"),
       height = 7, width = 11, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)


## Stats-------------------------------------------------
sap.risk.df <- readRDS("output/Cultivationrisk/df_total_sum_cache.rds") %>% 
  select(c(cell, cultrisk_en, count_percent)) %>% 
  pivot_wider(., names_from = cultrisk_en, values_from = count_percent) %>% 
  select(c(higher,cell))

# Mean
higher.mean <- sap.risk.df %>% summarise(mean(higher,na.rm=T)) %>% 
  pull()

# Forest coverage
risk.forestcover <- sap.risk %>% 
  project(., bay.f) %>%
  c(., bay.f) %>% #to make sure data is covering whole forest
  as.data.frame()
100-(sum(is.na(risk.forestcover$higher))/dim(risk.forestcover)[1]*100)


# Do the cells with low high cultivation risk have also low numbers of regeneration
quest <- readRDS("output/Cultivationrisk/df_total_sum_cache.rds") %>% 
  filter(cultrisk_en == "higher")
plot(log(quest$count_percent),log(quest$count_ha_cell))


# Which species drives the pattern?
df_total <- readRDS("output/Cultivationrisk/df_total_cache.rds")

df_total %>%
  group_by(species, cultrisk_en) %>%
  summarise(count_ha = sum(count_ha, na.rm = T), .groups = "drop") %>% # count_ha per species and cultrisk
  group_by(species) %>%
  mutate(count_ha_species = sum(count_ha),
         count_percent = count_ha/count_ha_species*100) %>% #percentage regeneration per cult risk category and species
  ungroup() %>%  
  filter(cultrisk_en == "higher") %>% # Percentage of higher risk!
  mutate(count_ha_higher = sum(count_ha),
         count_ha_higher_percent = (count_ha/count_ha_higher)*100) %>% 
  write.csv(., "output/Graphs/Tab_Cultivationrisk_driver.csv", row.names=F)


# BOX BAVARIA: TOTAL, DIVERSITY and CULT RISK -------------------------------------------------
## Cult risk ---------------------------------------------------------------
sap.risk <- rast("output/Cultivationrisk/Regeneration_cultivationrisk.tif") %>% 
  select(higher)
sap.risk.df.nona <- readRDS("output/Cultivationrisk/df_total_sum_cache.rds") %>% 
  select(c(cell, cultrisk_en, count_percent)) %>% 
  pivot_wider(., names_from = cultrisk_en, values_from = count_percent) %>% 
  select(c(higher,cell)) %>% 
  drop_na()

# Map
b.risk <-  
  ggplot()+
  theme_void()+
  geom_spatraster(data = sap.risk) +
  scale_fill_gradientn("Proportion of\nregeneration at\nhigh cultivation risk\n[%]",
                       colours = cult.col(50),
                       na.value = "transparent",
                       limits = c(0, 100)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = bavaria, fill = "transparent", colour = "black", linewidth = 0.1)

# Histogram
b.risk.h <-
  ggplot(sap.risk.df.nona, aes(x = higher)) +
  theme_classic() +
  geom_histogram(binwidth = 5,
                 boundary = 0,#-0.5
                 fill=cult.col(20), color="black") +
  xlab("Proportion of regeneration at\nhigh cultivation risk [%]") +
  ylab("Area [10\u00B3 ha]")+
  scale_y_continuous(limits=c(0,69e4),
                     breaks = c(0,2e5,4e5,6e5),
                     labels = c(0,200,400,600)) +
  theme(plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  annotate("segment", x = round(median(sap.risk.df.nona$higher, na.rm =T), digits=0), xend = 0, y = 65e4, yend = 65e4,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F") +
  annotate("text", x = round(median(sap.risk.df.nona$higher, na.rm =T))+ 2, y = 65e4, 
           hjust = 0, vjust = 0.5,
           label = paste(strwrap(paste0("50% of forest area has a cultivation risk <", round(median(sap.risk.df.nona$higher, na.rm =T), digits=1),"%"), 26), collapse = "\n"),
           colour = "#1D457F",
           size = 2) +
  annotate("segment", x = 100, xend = 75, y = 2e5, yend = 2e5,
           arrow = arrow(ends = "both", angle = 90, length = unit(.15,"cm")),
           colour = "#1D457F")+
  annotate("text", x = 100, y = 2.5e5, 
           hjust = 1, vjust = -0.1,
           label = paste(strwrap(paste0(round(table(sap.risk.df.nona$higher >= 75)[2]/dim(sap.risk.df.nona)[1]*100, digits=1),"% forest area has a high proportion of regeneration at risk (\u226575%)"), 30), collapse = "\n"),
           colour = "#1D457F",
           size = 2)

## Total Density ---------------------------------------------------------------------
regtot.bay <- rast("output/Graphs/Density.tif") %>% 
  project(.,sap.risk) %>% 
  crop(., sap.risk, mask = T)

# Categorize the raster values into classes
regtot.bay$class <- classify(regtot.bay, c(0, 1000, 2000, Inf), right= F)

b.tot <- 
  ggplot() + 
  theme_void() +
  geom_spatraster(data = regtot.bay$class)+
  scale_fill_manual(
    "Regeneration\ndensity [ha\u207B\u00B9]",
    na.translate = FALSE,
    labels = c("0-1000", "1000-2000", "\u22652000"),
    values = c("#9A3F07","#F2AF4AFF","grey90"),
    guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = bavaria, fill = "transparent", colour = "black", linewidth = 0.1)

regtot.bay.df <- as.data.frame(regtot.bay)
length(regtot.bay.df[regtot.bay.df$count_tot_ha <= 1000,]$count_tot_ha)/length(regtot.bay.df$count_tot_ha)*100


## Sp rich --------------------------------------------------------
sprich.bay <- rast("output/Graphs/Speciesrichness.tif") %>% 
  project(.,sap.risk) %>% 
  crop(., sap.risk, mask = T)
sprich.bay$class <- classify(sprich.bay$sprich, c(0, 2, 4, Inf))

b.sprich <- 
  ggplot() +
    theme_void() +
    geom_spatraster(data = sprich.bay$class) +
    scale_fill_manual(
      "Species richness",
      na.translate = FALSE,
      labels = c("0-2", "3-4", "\u22655"),
      values = c("#9A3F07","#F2AF4AFF","grey90"),
      guide = guide_legend(reverse = TRUE)) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 7)) +
  geom_spatvector(data = bavaria, fill = "transparent", colour = "black", linewidth = 0.1)

sprich.bay.df <- as.data.frame(sprich.bay)
mean(sprich.bay.df$sprich, na.rm=T)

# Stats
table(sprich.bay.df$sprich == 3 | sprich.bay.df$sprich == 4)[2]/length(sprich.bay.df$sprich)*100
table(sprich.bay.df$sprich >= 5)[2]/length(sprich.bay.df$sprich)*100

## All plots ---------------------------------------------------------------
b.tot + b.sprich + b.risk + b.risk.h +
  plot_layout(ncol = 2)+
  plot_annotation(tag_levels = "A")

ggsave(filename = paste0("output/Graphs/Regeneration_indicator_Bavaria.png"),
       height = 13, width = 17, units = "cm", dpi = 900,
       bg = "white",
       device=grDevices::png)

# TABLES ------------------------------------------------------------------
## Predictor variables overview---------------------------------------------
source("script/Model_vars.R")

Variables <-  Predictors::read_variables()
expl.vars <- c(fixed, random, spatial) %>%
  case_match(.,'northexp'~'northexp_der' # rename due to different data sources for modeling and predicting
             , 'eastexp'~ 'eastexp_der'
             , 'alt'~'alt_eudem'
             , 'blname'~'blname_der'
             , .default = .)

my.dat <- Variables %>%
  filter(varid %in% expl.vars)%>% 
  mutate(resolution = paste(resolution_x, "x", resolution_y, resolution_units)) %>% 
  select(c(domain, varid, description,unit_given, measurement_time, resolution,
           region, dsdir, varpath))

# Add new data
expl.vars[!(expl.vars %in% Predictors::availableVars())]

new.dat <- as.data.frame(
  matrix(ncol = dim(my.dat)[2], nrow=4,
         dimnames = list(NULL,names(my.dat)))) 

new.dat$varid <- c("wzp12_ba_ha_species", "yearmonth", "blname", "coordinates")
new.dat$unit_given <- c("m^2*ha^-1", "none", "none", "none")
new.dat$description <- c("basal area per ha of respective tree species",
                         "Month and Year of NFI measurement",
                         "federal state of Germany",
                         "coordinate of plot location")
new.dat$domain <- c("Stand structure", "Time","Space", "Space")
new.dat$region <- c("Germany", "Germany","Germany", "Germany")
new.dat$resolution <- c("none","none","none for NFI, 100 x 100 m for other data", "none")
new.dat$measurement_time <- c("2011-2012", "2011-2012", "2013 and 2022", "?")
new.dat

# combine data sets
my.dat <- rbind(my.dat, new.dat)

# save table
write.csv(my.dat,"output/Graphs/Explanatory_Variables_table.csv", row.names = F)


## Final species list ------------------------------------------------------
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv")

species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species

# All species that have good models n = 22
species.final. <- 
  species.tab %>% 
  filter(name.id %in% species.final)   

write.csv2(species.final., "output/Graphs/Regenerationdistribution_species.csv")


# All species that have good models and a cultivation risk map n = 17
species.final.cult <- 
  species.tab %>% 
  filter(name.id %in% species.final) %>% 
  filter(!is.na(name.cultrisk)) 

write.csv2(species.final.cult, "output/Graphs/Cultivationrisk_species.csv")


## Model and density summary -----------------------------------------------------------
# Model performance
model.all <- read.csv2("output/Fits/Sapling/h50d7_Germany/Sapling_model_summary.csv") %>% 
  select(c(species, rsq.test.median, mae.relative.median)) %>% 
  mutate(mae.rule = ifelse(mae.relative.median <= 2, "1", "0"),
         rsq.rule = ifelse(rsq.test.median >= 0.1, "1", "0"),
         perfcrit = ifelse(mae.rule ==1 & rsq.rule ==1, "1","0"))

# Total  map means
mean.maps <- read.csv2("output/Graphs/Regeneration_density_stats.csv") %>% 
  rename(name.id = X)

model.merge <- merge(model.all, mean.maps, by.x = "species", by.y = "name.id", all = T)

# Availability of cult risk maps
species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv") %>% 
  mutate(cultavail = ifelse(is.na(name.cultrisk), 0,1)) %>% 
 select(c(name.id, name.scient, cultavail)) 

model.merge <- merge(model.merge, species.tab, by.x = "species", by.y = "name.id", all.x = T)

write.csv(model.merge, "output/Graphs/Tab_ModelMap_summary.csv",row.names = F)

## Map covering forest area----------------------------

### Germany -----------------------------------------------------------------
dt.f <- rast("data/ForestArea_germany_100m.tif")
map <-  rast(paste0("output/Predictions/Regeneration_Abies.alba.tif"))

cover.dt <- c(dt.f,map$Abies.alba) %>% 
  as.data.frame()

table(is.na(cover.dt$Abies.alba))/dim(cover.dt)[1]*100


### Bavaria -----------------------------------------------------------------
bay.f <- rast("data/ForestArea_bavaria_100m.tif")
map <-  rast(paste0("output/Predictions/Regeneration_Abies.alba.tif"))
crs(bay.f)==crs(map)

cover.bay <- map$Abies.alba %>% 
  crop(., bay.f, mask=T) %>%
  c(., bay.f) %>% #to make sure data is covering whole forest
  as.data.frame()
table(is.na(cover.bay$Abies.alba))/dim(cover.bay)[1]*100

## BWI regeneration stats  ------------------------------------------------------
### Data prep --------------------------------------------------------------------
# Species 
source("script/Species_select.R") # all species species.vect

species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species

species.final.cult <- read.csv("data/DE_BWI3_regeneration_species.csv") %>% 
  filter(name.id %in% species.final) %>% 
  filter(!is.na(name.cultrisk)) %>% 
  select(name.id) %>% 
  unlist(as.vector(.))

# Regeneration
reg <- readRDS("data/DE_BWI3_regeneration_h50d7.rds") %>%
  unnest(data) %>%
  mutate(tax= droplevels(tax)) %>% 
  select(c(tax,plotid, count, countarea))

bl <- readRDS("data/DE_BWI3_explvars_df.rds") %>%
  select(c(plotid, blname_loc))

# join Bundesland Information
reg <-  inner_join(reg,bl, by = "plotid")

# filter Bavaria
reg.by <- reg %>% 
  filter(blname_loc == "BY") 

# Output data frame
df <- data.frame(matrix(ncol = 4, nrow = 2))
colnames(df) <- c("count_tot", "count_selsp","count_validmod","count_validmod_cult")
rownames(df) <- c("Germany", "Bavaria")


### How much mean regeneration per ha ---------------------------------------------------------
BWI_reg <- readRDS("data/DE_BWI3_regeneration_h50d7.rds")

BWI_regtot <- 
  BWI_reg %>%   
  unnest(data) %>% 
  select(c(tax,count,plotid)) %>% 
  as.data.frame() %>% 
  group_by(plotid) %>% 
  summarise(count_tot = sum(count)) %>% 
  mutate(count_tot_ha = (count_tot*10000)/(2^2*pi))

E <- readRDS("data/DE_BWI3_explvars_df.rds") 

merg <- merge(BWI_regtot, E)

dat = merg %>% filter(blname_loc=="BY")
summary(dat$count_tot)
sd(dat$count_tot)


### How much BWI regeneration count data is covered by our study?  ---------------------------------
# by selected species
# by species of valid models
# by species of valid models and available cultivation risk data

#Total count
df["Germany","count_tot"] <- reg %>%
  summarise(sum(count)) 

df["Bavaria","count_tot"] <- reg.by %>%
  summarise(sum(count)) 

# Selected species 
df["Germany","count_selsp"] <- reg %>%
  filter(tax %in% species.vect) %>% 
  summarise(sum(count))

df["Bavaria","count_selsp"] <- reg.by %>%
  filter(tax %in% species.vect) %>% 
  summarise(sum(count)) 

# Validated model 
df["Germany","count_validmod"] <- reg %>%
  filter(tax %in% species.final) %>% 
  summarise(sum(count))

df["Bavaria","count_validmod"] <- reg.by %>%
  filter(tax %in% species.final) %>% 
  summarise(sum(count)) 

# Cultivation risk 
df["Bavaria","count_validmod_cult"] <- reg.by %>%
  filter(tax %in% species.final.cult) %>% 
  summarise(sum(count)) 

# Percentages 
df %<>% 
  mutate(count_validmod_per = (count_validmod/count_selsp)*100,
         count_validmod_cult_per = (count_validmod_cult/count_selsp)*100)
df


# SUPPORTING INFORMATION----------------------------------------------------------------

## Model spat autocorr and blocked cross validation ------------------------------------------------
model.autocorr.cv <- read.csv2("output/Fits/Sapling/h50d7_Germany/Sapling_model_summary.csv") %>% 
  select(c(species, spatautocorr.Iobserved, spatautocorr.p,cv.sp.range, cv.set.range, cv.blocknr))

species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv") %>% 
  select(c(name.id, name.scient)) 

merge(model.autocorr.cv, species.tab, by.x = "species", by.y = "name.id", all.x = T) %>% 
  write.csv(., "output/Graphs/Tab_Model_spatcorr_cv.csv", row.names=F)


## Mean species regenration densities BWI? --------------------------------
BWI_reg <- readRDS("data/DE_BWI3_regeneration_h50d7.rds")
source("script/Species_select.R")

species.sel = as.character(species.vect) # in case "Sorbus.domestica" is not part of it
ifelse("Sorbus.domestica" %in% species.sel, species.sel <- species.sel, species.sel <- c(species.sel,"Sorbus.domestica"))

species.tab <- read.csv("data/DE_BWI3_regeneration_species.csv") %>% 
  select(c(name.id, name.scient)) 
species.final <- readRDS("output/Fits/Sapling/h50d7_Germany/Sapling_model_final.rds")$species

BWI_specdens <- BWI_reg %>% 
  unnest(data) %>% # turn multilayer tibble into long data.frame
  filter(tax %in% species.sel) %>% # select only tree species we are interested in
  mutate(count_ha = (count*10000)/(2^2*pi), # convert to ha
         tax = droplevels(tax)) %>% # drop levels with no occurrence
  select(c(tax, count_ha, plotid)) %>% # select only relevant columns
  group_by(tax) %>% # for each species
  summarise(sum_count_ha = sum(count_ha), # total regeneration sum
            mean_count_ha = mean(count_ha), # mean density
            sd_count_ha = sd(count_ha)) %>% # sd denstiy
  ungroup() %>% 
  mutate(tot_perc_count_ha = sum_count_ha/sum(sum_count_ha)*100) # proportion of total regeneration

BWI_specdens %>% 
  mutate(densmap = ifelse(tax %in% species.final, 1, 0)) %>% # if density map exists 1, if not 0, mean model performance criteria met
  merge(., species.tab, by.x = "tax", by.y = "name.id", all.x = T) %>% # add nice scientific species names
  write.csv(., "output/Graphs/Tab_BWI_reg_summary.csv", row.names=F)

